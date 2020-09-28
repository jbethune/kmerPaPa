use anyhow::Result;
use clap::{App, Arg};

use pattern_partition_prediction::PaPaPred;
use twobit::TwoBitFile;

mod compare;
mod counts;
mod enumerate;
mod error;
mod expect;
mod io;
mod observed;
mod sample;
mod transform;

type Float = f32;
type MutationType = mutexpect::MutationType;

use crate::compare::compare_mutations;
use crate::enumerate::enumerate_possible_mutations;
use crate::error::MissingCommandLineArgumentError;
use crate::expect::expected_number_of_mutations;
use crate::observed::classify_mutations;
use crate::observed::read_mutations_from_file as read_observed_mutations_from_file;
use crate::sample::sample_mutations;

fn require_initialization<'a, T>(
    value: &'a Option<T>,
    cli_argument: &'static str,
) -> std::result::Result<&'a T, MissingCommandLineArgumentError> {
    match value {
        Some(v) => Ok(v),
        None => Err(MissingCommandLineArgumentError::new(cli_argument)),
    }
}

fn main() -> Result<()> {
    let app = App::new("genovo")
        .version("0.1.3")
        .author("Jörn Bethune")
        .about("Determine genes enriched with de-novo mutations")
        .after_help("If no --action is given, all actions are executed.\n\
                     Possible actions are: transform, enumerate, expect, sample, classify, compare" )
        .arg(Arg::with_name("action")
             .long("action")
             .value_name("ACTION")
             .help("Only run a specific step in the pipeline")
             .takes_value(true))

        //raw data arguments
        .arg(Arg::with_name("gff3")
             .long("gff3")
             .value_name("FILE")
             .help("gff3 gene annotations file")
             .takes_value(true))
        .arg(Arg::with_name("genome")
             .long("genome")
             .value_name("FILE")
             .help("A 2bit reference genome sequence file")
             .takes_value(true))
        .arg(Arg::with_name("mutation-probabilities")
             .long("mutation-probabilities")
             .value_name("FILE")
             .help("A pattern partition prediction mutation probability table")
             .takes_value(true))
        .arg(Arg::with_name("observed-mutations")
             .long("observed-mutations")
             .value_name("FILE")
             .help("A vcf-like file containing observed point mutations")
             .takes_value(true))

        // input/output file arguments
        .arg(Arg::with_name("genomic-regions")
             .long("genomic-regions")
             .value_name("FILE")
             .help("Locations of exons, CDS and their phases for each gene")
             .takes_value(true))
        .arg(Arg::with_name("possible-mutations")
             .long("possible-mutations")
             .value_name("FILE")
             .help("A list of all possible point mutations for each gene")
             .takes_value(true))
        .arg(Arg::with_name("classified-mutations")
             .long("classified-mutations")
             .value_name("FILE")
             .help("Observed, classified point mutations")
             .takes_value(true))
        .arg(Arg::with_name("expected-mutations")
             .long("expected-mutations")
             .value_name("FILE")
             .help("Expected number of point mutations per gene")
             .takes_value(true))
        .arg(Arg::with_name("sampled-mutations")
             .long("sampled-mutations")
             .value_name("FILE")
             .help("Sampled number of point mutations per gene")
             .takes_value(true))
        .arg(Arg::with_name("significant-mutations")
             .long("significant-mutations")
             .value_name("FILE")
             .help("Statistical test results for every gene")
             .default_value("-")
             .takes_value(true))

        // non-file args
        .arg(Arg::with_name("id")
             .long("--id")
             .value_name("NAME")
             .help("Only process a gene/transcript with the given ID (useful for parallel processing)")
             .takes_value(true))
        .arg(Arg::with_name("number-of-random-samples")
             .long("--number-of-random-samples")
             .value_name("NUMBER")
             .default_value("1000")
             .help("The number of random samples that should be generated")
             .takes_value(true))
    ;

    let matches = app.get_matches();
    let run_all = matches.value_of("action").is_none();
    let id = matches.value_of("id");

    /*
     * Depending on what --action is given, each step may or may not produce results.
     * Therefore the variables are all Option's.
     */

    let ref_genome = {
        if let Some(ref_genome_file) = matches.value_of("genome") {
            Some(TwoBitFile::open(ref_genome_file, false)?)
        } else {
            None
        }
    };

    let papa = {
        if let Some(papa_file) = matches.value_of("mutation-probabilities") {
            Some(PaPaPred::new(papa_file, Some(5))?) // 5 to have 2 flanking bases around a point mutation
        } else {
            None
        }
    };

    let observed_mutations = {
        if let Some(observed_mutations_file) = matches.value_of("observed-mutations") {
            Some(read_observed_mutations_from_file(
                observed_mutations_file,
                -1,
            )?) //TODO expose adjustment parameter to CLI
        } else {
            None
        }
    };

    // action=transform
    let regions = {
        if run_all || matches.value_of("action") == Some("transform") {
            if let Some(gff3) = matches.value_of("gff3") {
                let regions = transform::transform_gff3_annotations(gff3, id)?;
                if let Some(regions_file) = matches.value_of("genomic-regions") {
                    transform::write_sequence_annotations_to_file(regions_file, &regions)?;
                }
                if !run_all {
                    // we are done here
                    return Ok(());
                }
                Some(regions)
            } else {
                return Err(anyhow::anyhow!("Please provide the --gff3 parameter"));
            }
        } else if let Some(regions_file) = matches.value_of("genomic-regions") {
            Some(mutexpect::read_sequence_annotations_from_file(
                regions_file,
                id,
            )?)
        } else {
            None
        }
    };

    //action=enumerate
    let possible_mutations = {
        if run_all || matches.value_of("action") == Some("enumerate") {
            let possible_mutations = enumerate_possible_mutations(
                require_initialization(&regions, "--genomic-regions")?,
                require_initialization(&ref_genome, "--genome")?,
                require_initialization(&papa, "--mutation-probabilities")?,
                true,
                id,
            )?;

            if let Some(possible_mutations_file) = matches.value_of("possible-mutations") {
                enumerate::write_to_file(possible_mutations_file, &possible_mutations)?
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(possible_mutations)
        } else if let Some(possible_mutations_file) = matches.value_of("possible-mutations") {
            Some(enumerate::read_from_file(possible_mutations_file)?)
        } else {
            None
        }
    };

    let expected_mutations = {
        if run_all || matches.value_of("action") == Some("expect") {
            let expected_mutations = expected_number_of_mutations(
                require_initialization(&possible_mutations, "--possible-mutations")?,
                id,
            )?;
            if let Some(expected_mutations_file) = matches.value_of("expected-mutations") {
                expect::write_to_file(expected_mutations_file, &expected_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(expected_mutations)
        } else if let Some(expected_mutations_file) = matches.value_of("expected-mutations") {
            Some(expect::read_from_file(expected_mutations_file)?)
        } else {
            None
        }
    };

    let sampled_mutations = {
        if run_all || matches.value_of("action") == Some("sample") {
            let iterations: usize = matches
                .value_of("number-of-random-samples")
                .expect("clap default value")
                .parse()
                .unwrap(); //TODO proper error handling
            let sampled_mutations = sample_mutations(
                require_initialization(&possible_mutations, "--possible-mutations")?,
                iterations,
                true,
                id,
            )?;

            if let Some(sampled_mutations_file) = matches.value_of("sampled-mutations") {
                sample::write_to_file(sampled_mutations_file, &sampled_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(sampled_mutations)
        } else if let Some(sampled_mutations_file) = matches.value_of("sampled-mutations") {
            Some(sample::read_from_file(sampled_mutations_file)?)
        } else {
            None
        }
    };

    std::mem::drop(possible_mutations); // let's free up some memory

    let classified_mutations = {
        if run_all || matches.value_of("action") == Some("classify") {
            let classified_mutations = classify_mutations(
                require_initialization(&observed_mutations, "--observed-mutations")?,
                require_initialization(&regions, "--genomic-regions")?,
                require_initialization(&ref_genome, "--genome")?,
                id,
            )?;

            if let Some(classified_mutations_file) = matches.value_of("classified-mutations") {
                observed::write_to_file(classified_mutations_file, &classified_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(classified_mutations)
        } else if let Some(classified_mutations_file) = matches.value_of("classified-mutations") {
            Some(observed::read_from_file(classified_mutations_file)?)
        } else {
            None
        }
    };

    let _significant_mutations = {
        if run_all || matches.value_of("action") == Some("compare") {
            let significant_mutations = compare_mutations(
                require_initialization(&classified_mutations, "--classified-mutations")?,
                require_initialization(&expected_mutations, "--expected-mutations")?,
                require_initialization(&sampled_mutations, "--sampled-mutations")?,
                id,
            )?;

            if let Some(significant_mutations_file) = matches.value_of("significant-mutations") {
                compare::write_to_file(significant_mutations_file, &significant_mutations)?;
            }
            if !run_all {
                // we are done here
                return Ok(());
            }
            Some(significant_mutations)
        } else {
            // last step in the pipeline. No work needs to be done
            None
        }
    };

    if !run_all {
        return Err(anyhow::anyhow!(
            "Invalid --action parameter: {}",
            matches.value_of("action").expect("clap")
        ));
    }
    Ok(())
}
