//! This module contains the FunctionAggregator struct that is responsible for aggregating the
//! functional annotations of proteins.

use std::collections::{HashMap, HashSet};
use serde::Serialize;


use crate::proteins::Protein;

/// A struct that represents the functional annotations once aggregated
#[derive(Debug, Serialize)]
pub struct FunctionalAggregation {
    /// A HashMap representing how many GO, EC and IPR terms were found
    pub counts: HashMap<String, usize>,
    /// A HashMap representing how often a certain functional annotation was found
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
    pub fn aggregate(&self, proteins: Vec<&Protein>) -> FunctionalAggregation {
        // Keep track of the proteins that have a certain annotation
        let mut proteins_with_ec: HashSet<String> = HashSet::new();
        let mut proteins_with_go: HashSet<String> = HashSet::new();
        let mut proteins_with_ipr: HashSet<String> = HashSet::new();

        // Keep track of the counts of the different annotations
        let mut data: HashMap<String, u32> = HashMap::new();

        for protein in proteins.iter() {
            for annotation in protein.get_functional_annotations().split(';') {
                match annotation.chars().next() {
                    Some('E') => proteins_with_ec.insert(protein.uniprot_id.clone()),
                    Some('G') => proteins_with_go.insert(protein.uniprot_id.clone()),
                    Some('I') => proteins_with_ipr.insert(protein.uniprot_id.clone()),
                    _ => false
                };

                data.entry(annotation.to_string()).and_modify(|c| *c += 1).or_insert(1);
            }
        }

        let mut counts: HashMap<String, usize> = HashMap::new();
        counts.insert("all".to_string(), proteins.len());
        counts.insert("EC".to_string(), proteins_with_ec.len());
        counts.insert("GO".to_string(), proteins_with_go.len());
        counts.insert("IPR".to_string(), proteins_with_ipr.len());

        data.remove("");

        FunctionalAggregation { counts, data }
    }


    /// Aggregates the functional annotations of proteins
    ///
    /// # Arguments
    /// * `proteins` - A vector of proteins
    ///
    /// # Returns
    ///
    /// Returns a list of lists with all the functional annotations per protein
    pub fn get_all_functional_annotations(&self, proteins: &[&Protein]) -> Vec<Vec<String>> {
        proteins
            .iter()
            .map(
                |&prot| 
                    prot.get_functional_annotations()
                        .split(';')
                        .map(|ann| ann.to_string())
                        .filter(|s| !s.is_empty())
                        .collect()
            )
            .collect::<Vec<Vec<String>>>()
    }
}
