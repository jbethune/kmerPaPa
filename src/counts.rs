use std::convert::TryFrom;
use std::iter::{IntoIterator, Iterator};
use std::ops::{Add, AddAssign};

use crate::{Float, MutationType};

#[derive(Copy, Clone, Default, Debug, PartialEq)]
pub struct MutationTypeCounts<Count> {
    pub unknown: Count,
    pub synonymous: Count,
    pub missense: Count,
    pub nonsense: Count,
    pub start_codon: Count,
    pub stop_loss: Count,
    pub splice_site: Count,
}

impl<Count: Copy + Add<Output = Count> + AddAssign> MutationTypeCounts<Count> {
    pub fn get(&self, mutation_type: MutationType) -> Count {
        match mutation_type {
            MutationType::Unknown => self.unknown,
            MutationType::Synonymous => self.synonymous,
            MutationType::Missense => self.missense,
            MutationType::Nonsense => self.nonsense,
            MutationType::StartCodon => self.start_codon,
            MutationType::StopLoss => self.stop_loss,
            MutationType::SpliceSite => self.splice_site,
        }
    }

    pub fn add(&mut self, mutation_type: MutationType, value: Count) {
        match mutation_type {
            MutationType::Unknown => self.unknown += value,
            MutationType::Synonymous => self.synonymous += value,
            MutationType::Missense => self.missense += value,
            MutationType::Nonsense => self.nonsense += value,
            MutationType::StartCodon => self.start_codon += value,
            MutationType::StopLoss => self.stop_loss += value,
            MutationType::SpliceSite => self.splice_site += value,
        }
    }

    pub fn mutation_types() -> Vec<MutationType> {
        vec![
            MutationType::Unknown,
            MutationType::Synonymous,
            MutationType::Missense,
            MutationType::Nonsense,
            MutationType::StartCodon,
            MutationType::StopLoss,
            MutationType::SpliceSite,
        ]
    }
}

pub type ExpectedMutationCounts = MutationTypeCounts<Float>;
pub type ObservedMutationCounts = MutationTypeCounts<usize>;

impl<Count: Clone> IntoIterator for MutationTypeCounts<Count> {
    type Item = Count;
    type IntoIter = MutationTypeCountsIter<Count>;
    fn into_iter(self) -> Self::IntoIter {
        Self::IntoIter::new(self)
    }
}

pub struct MutationTypeCountsIter<Count> {
    counts: MutationTypeCounts<Count>,
    index: usize,
}

impl<Count> MutationTypeCountsIter<Count> {
    fn new(counts: MutationTypeCounts<Count>) -> Self {
        Self { counts, index: 0 }
    }
}

impl<Count: Clone> Iterator for MutationTypeCountsIter<Count> {
    type Item = Count;

    fn next(&mut self) -> Option<Count> {
        let result = match self.index {
            0 => &self.counts.unknown,
            1 => &self.counts.synonymous,
            2 => &self.counts.missense,
            3 => &self.counts.nonsense,
            4 => &self.counts.start_codon,
            5 => &self.counts.stop_loss,
            6 => &self.counts.splice_site,
            _ => return None,
        }
        .clone();
        self.index += 1;
        Some(result)
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct DefaultCounter {
    values: Vec<usize>,
}

impl DefaultCounter {
    pub fn new() -> DefaultCounter {
        DefaultCounter { values: Vec::new() }
    }

    pub fn inc(&mut self, index: usize) {
        while self.values.len() <= index {
            self.values.push(0)
        }
        self.values[index] += 1
    }

    /// Calculate the p-value for observed a given number of events (represented by the vector
    /// index).
    ///
    /// This is a one-sided test with the right tail being significant
    pub fn p_values(&self) -> PValues {
        let mut result = Vec::new();
        let total: Float = self.values.iter().sum::<usize>() as Float;
        let mut accumulator: Float = 0.0;

        for count in self.values.iter().rev() {
            // I go from right to left, because it's easier to calculate
            accumulator += *count as Float;
            result.push(accumulator / total);
        }

        result.reverse();
        PValues { p_values: result }
    }
}

impl ToString for DefaultCounter {
    fn to_string(&self) -> String {
        let mut result = String::with_capacity(self.values.len() * 4); // 4 is a guestimate
        for value in &self.values {
            if !result.is_empty() {
                result.push('|');
            }
            result.push_str(&value.to_string());
        }
        result
    }
}

impl TryFrom<&str> for DefaultCounter {
    type Error = std::num::ParseIntError;
    fn try_from(string: &str) -> std::result::Result<Self, Self::Error> {
        let mut result = Self::new();
        for str_count in string.split('|') {
            result.values.push(str_count.parse()?);
        }
        Ok(result)
    }
}

pub struct PValues {
    p_values: Vec<Float>,
}

impl PValues {
    pub fn n_hits_or_more(&self, n: usize) -> Float {
        self.p_values.get(n).copied().unwrap_or(0.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_p_values() {
        let mut counter = DefaultCounter::new();
        counter.inc(5);
        for _ in 0..10 {
            counter.inc(6);
        }
        for _ in 0..20 {
            counter.inc(7);
        }
        for _ in 0..30 {
            counter.inc(8);
        }
        for _ in 0..39 {
            counter.inc(9);
        }
        let p_values = counter.p_values();
        assert_eq!(p_values.n_hits_or_more(0), 1.0);
        assert_eq!(p_values.n_hits_or_more(5), 1.0);
        assert_eq!(p_values.n_hits_or_more(6), 99.0 / 100.0);
        assert_eq!(p_values.n_hits_or_more(7), 89.0 / 100.0);
        assert_eq!(p_values.n_hits_or_more(8), 69.0 / 100.0);
        assert_eq!(p_values.n_hits_or_more(9), 39.0 / 100.0);
        assert_eq!(p_values.n_hits_or_more(10), 0.0);
    }
}
