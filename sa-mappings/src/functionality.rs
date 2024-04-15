//! This module contains the FunctionAggregator struct that is responsible for aggregating the
//! functional annotations of proteins.

use std::collections::{HashMap, HashSet};

use serde_json::json;

use crate::proteins::Protein;

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
    pub fn aggregate(&self, proteins: Vec<String>) -> String {
        let mut data: HashMap<String, u32> = HashMap::new();

        let mut ec_numbers: HashSet<String> = HashSet::new();
        let mut go_terms: HashSet<String> = HashSet::new();
        let mut interpros: HashSet<String> = HashSet::new();

        for fa in proteins.iter() {
            for annotation in fa.split(';') {
                let annotation_string = annotation.to_string();

                match annotation_string.chars().next() {
                    Some('E') => ec_numbers.insert(annotation_string.clone()),
                    Some('G') => go_terms.insert(annotation_string.clone()),
                    Some('I') => interpros.insert(annotation_string.clone()),
                    _ => false
                };

                data.entry(annotation_string).and_modify(|e| *e += 1).or_insert(1);
            }
        }

        json!({
            "num": {
                "all": proteins.len(),
                "EC": ec_numbers.len(),
                "GO": go_terms.len(),
                "IPR": interpros.len(),
            },
            "data": data
        }).to_string()
    }
}
