use std::io::{BufWriter, Write};
use std::path::Path;

use anyhow::{Context, Result};
use serde::{Deserialize, Serialize};

use mutexpect::interval::Interval;
use mutexpect::{MutationType, PointMutationClassifier};
use tabfile::Tabfile;
use twobit::TwoBitFile;

use crate::compare::tally_up_observed_mutations;
use crate::counts::ObservedMutationCounts;
use crate::error::ParseError;
use crate::io::{get_reader, get_writer};

pub fn classify_mutations(
    observed_mutations: &[PointMutation],
    annotations: &[mutexpect::SeqAnnotation],
    genome: &TwoBitFile,
    filter_for_id: Option<&str>,
) -> Result<Vec<AnnotatedPointMutation>> {
    let mut result = Vec::new();

    let flank = 2; // number of flanking bases left and right needed to classify all coding point mutations

    for annotation in annotations {
        if let Some(id) = filter_for_id {
            if id != annotation.name {
                continue;
            }
        }

        let seq_of_region: Vec<char> = genome
            .sequence(
                &annotation.chr,
                annotation.range.start - flank,
                annotation.range.stop + flank,
            )?
            .chars()
            .collect();
        assert_eq!(seq_of_region.len(), 2 * flank + annotation.range.len());

        let classifier = PointMutationClassifier::new(&annotation, 2);
        let mut relevant_mutations =
            filter_observed_mutations(&observed_mutations, &annotation.chr, annotation.range);
        for mutation in &mut relevant_mutations {
            let sequence_context: Vec<char> = {
                assert!(annotation.range.start <= mutation.position);
                let middle = mutation.position - annotation.range.start + flank;
                seq_of_region[middle - flank..middle + flank + 1].into()
            };
            assert_eq!(sequence_context[2], mutation.from); // sanity-check right reference genome

            let overlapping_intron = annotation.find_intron(mutation.position);

            let mut mutation_type = classifier.classify_by_position(
                mutation.position,
                &sequence_context,
                &overlapping_intron, // may be none
            );

            if mutation_type == MutationType::Unknown {
                if let Some(overlapping_cds) = annotation.find_cds(mutation.position) {
                    mutation_type = classifier.classify_coding_mutation(
                        mutation.position,
                        &sequence_context,
                        mutation.to,
                        &overlapping_cds,
                    );
                }
            }
            result.push(mutation.annotate(annotation.name.clone(), mutation_type));
        }
    }
    Ok(result)
}

pub fn read_mutations_from_file<P: AsRef<Path>>(
    file: P,
    adjust: i64,
) -> Result<Vec<PointMutation>> {
    let mut result = Vec::new();

    let tsv_reader = Tabfile::open(&file)
        .with_context(|| format!("failed to open file {}", &file.as_ref().display()))?
        .comment_character('#')
        .separator(' ');

    for record_result in tsv_reader {
        let record = record_result?;
        let fields = record.fields();
        if fields.len() < 4 {
            return Err( ParseError::new(format!("Bad format in line {}. Expecting at least 4 tab-delimited fields: chr, pos, ref, alt", record.line_number()) ).into());
        }
        if fields[2].len() != 1 || fields[3].len() != 1 {
            continue; // skip indels
        }
        let chromosome = fields[0].to_string();
        let position = {
            let value = fields[1].parse::<i64>()?;
            (value + adjust) as usize
        };
        let from = fields[2].chars().next().unwrap();
        let to = fields[3].chars().next().unwrap();
        result.push(PointMutation {
            chromosome,
            position,
            from,
            to,
        });
    }

    Ok(result)
}

pub struct PointMutation {
    chromosome: String,
    position: usize,
    from: char,
    to: char,
}

impl PointMutation {
    pub fn annotate(
        &self,
        region_name: String,
        mutation_type: MutationType,
    ) -> AnnotatedPointMutation {
        AnnotatedPointMutation {
            region_name,
            mutation_type,
            chromosome: self.chromosome.clone(),
            position: self.position,
            from: self.from,
            to: self.to,
        }
    }
}

#[derive(Debug, Serialize, Deserialize, PartialEq, Eq)]
pub struct AnnotatedPointMutation {
    pub region_name: String,
    pub mutation_type: MutationType,
    pub chromosome: String,
    pub position: usize,
    pub from: char,
    pub to: char,
}

fn filter_observed_mutations<'a>(
    mutations: &'a [PointMutation],
    chr: &str,
    genomic_region: Interval,
) -> Vec<&'a PointMutation> {
    let mut result = Vec::new();
    for mutation in mutations {
        // I assume no particular ordering. Otherwise a binary search might be faster
        if mutation.chromosome == chr && genomic_region.contains(mutation.position) {
            result.push(mutation)
        }
    }
    result
}

// serialization stuff //

pub fn write_to_file(out_path: &str, annotated_mutations: &[AnnotatedPointMutation]) -> Result<()> {
    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut csv_writer = csv::WriterBuilder::new()
        .delimiter(b'\t')
        .from_writer(writer);
    for mutation in annotated_mutations {
        csv_writer.serialize(mutation)?;
    }
    Ok(())
}

/// Write a file with the following format:
/// ```
/// transcript_id<tab>synonymous<tab>missense<tab>...
/// ENSG1234<tab>synonymous_mutation_count<tab>missense_mutation_count<tab>...
/// ```
///
pub fn sum_up_and_write_to_file(
    out_path: &str,
    annotated_mutations: &[AnnotatedPointMutation],
) -> Result<()> {
    let transcript_mutation_counts = tally_up_observed_mutations(annotated_mutations, None);

    let writer = get_writer(out_path)
        .with_context(|| format!("failed to open file {} for writing", out_path))?;
    let mut buf_writer = BufWriter::new(writer);
    buf_writer.write_all(b"name")?;
    for mut_type in ObservedMutationCounts::mutation_types() {
        buf_writer.write_fmt(format_args!("\t{}", mut_type))?;
    }
    buf_writer.write_all(b"\n")?; // end of header line

    // for every transcript
    for (name, counts) in transcript_mutation_counts {
        buf_writer.write(name.as_bytes())?;
        for count in counts {
            buf_writer.write_fmt(format_args!("\t{}", count))?;
        }
        buf_writer.write_all(b"\n")?;
    }
    Ok(())
}

pub fn read_from_file(in_path: &str) -> Result<Vec<AnnotatedPointMutation>> {
    let mut result = Vec::new();
    let reader = get_reader(in_path)
        .with_context(|| format!("failed to open file {} for reading", in_path))?;
    let mut csv_reader = csv::ReaderBuilder::new()
        .delimiter(b'\t')
        .from_reader(reader);
    for row_result in csv_reader.deserialize() {
        let mutation: AnnotatedPointMutation = row_result?;
        result.push(mutation);
    }
    Ok(result)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_observed_mutations_io() {
        let path = "/tmp/unit_test.observed_mutations";
        let mut om = Vec::new();
        write_to_file(path, &om).unwrap();
        let om2 = read_from_file(path).unwrap();
        assert_eq!(om2.len(), 0);

        om.push(AnnotatedPointMutation {
            region_name: "my_gene".to_string(),
            mutation_type: MutationType::Synonymous,
            chromosome: "chr42".to_string(),
            position: 42,
            from: 'A',
            to: 'T',
        });
        write_to_file(path, &om).unwrap();
        let om2 = read_from_file(path).unwrap();
        assert_eq!(om, om2);

        om.push(AnnotatedPointMutation {
            region_name: "my_other_gene".to_string(),
            mutation_type: MutationType::Nonsense,
            chromosome: "chrM".to_string(),
            position: 4,
            from: 'C',
            to: 'G',
        });

        write_to_file(path, &om).unwrap();
        let om2 = read_from_file(path).unwrap();
        assert_eq!(om, om2);
    }
}
