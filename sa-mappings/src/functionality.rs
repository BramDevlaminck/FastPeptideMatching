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

        for fa in proteins.iter() {
            for annotation in fa.split(';') {
                let annotation = annotation.to_string();

                // Update the counts
                match annotation.chars().next() {
                    Some('E') => counts_array[0] += 1,
                    Some('G') => counts_array[1] += 1,
                    Some('I') => counts_array[2] += 1,
                    _ => ()
                };

                data.entry(annotation).and_modify(|e| *e += 1).or_insert(1);
            }
        }

        counts.insert("all".to_string(), counts_array.iter().sum());
        counts.insert("EC".to_string(), counts_array[0]);
        counts.insert("GO".to_string(), counts_array[1]);
        counts.insert("IPR".to_string(), counts_array[2]);

        FunctionalAggregation { counts, data }
    }
}
