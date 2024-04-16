//! This module contains the FunctionAggregator struct that is responsible for aggregating the
//! functional annotations of proteins.

use std::collections::{HashMap, HashSet};

use serde::Serialize;
use serde_json::json;

use crate::proteins::Protein;

#[derive(Debug, Serialize)]
pub struct FunctionalAggregation {
    pub counts: HashMap<String, u32>,
    pub data: HashMap<String, u32>,
}

/// A struct that represents a function aggregator
pub struct FunctionAggregator {}

impl FunctionAggregator {
    /// Aggregates the functional annotations of proteins
    ///
    /// # Arguments
    /// * `proteins` - A vector of proteins
    ///
    /// # Returns
    ///
    /// Returns a JSON string containing the aggregated functional annotations
    pub fn aggregate(&self, proteins: Vec<String>) -> FunctionalAggregation {
        let mut counts_array: [u32; 3] = [0; 3];

        let mut counts: HashMap<String, u32> = HashMap::new();
        let mut data: HashMap<String, u32> = HashMap::new();

        let mut update = |i: usize, fa: &str| {
            counts_array[i] += 1;
            data.entry(fa.to_string()).and_modify(|e| *e += 1).or_insert(1);
        };

        for fa in proteins.iter() {
            for annotation in fa.split(';') {
                match annotation.chars().next() {
                    Some('E') => update(0, annotation),
                    Some('G') => update(1, annotation),
                    Some('I') => update(2, annotation),
                    _ => ()
                };
            }
        }

        counts.insert("all".to_string(), counts_array.iter().sum());
        counts.insert("EC".to_string(), counts_array[0]);
        counts.insert("GO".to_string(), counts_array[1]);
        counts.insert("IPR".to_string(), counts_array[2]);

        FunctionalAggregation { counts, data }
    }
}
