use serde::Serialize;
use std::collections::HashMap;
use tsv_utils::Protein;

#[derive(Debug, Serialize)]
pub struct FunctionalAnnotations {
    ec: HashMap<String, u64>,
    go: HashMap<String, u64>,
    interpro: HashMap<String, u64>,
}

impl FunctionalAnnotations {
    pub fn get_uniprot_ids(proteins: &[&Protein]) -> Vec<String> {
        proteins
            .iter()
            .map(|&protein| protein.uniprot_id.clone())
            .collect()
    }

    fn count_occurrences(mapping: &mut HashMap<String, u64>, values: &Vec<String>) {
        for val in values {
            mapping
                .entry(val.clone())
                .and_modify(|v| *v += 1)
                .or_insert(1);
        }
    }

    pub fn new(proteins: &[&Protein]) -> Self {
        let mut ec_aggregate = HashMap::new();
        let mut go_aggregate = HashMap::new();
        let mut interpro_aggregate = HashMap::new();

        for &protein in proteins {
            Self::count_occurrences(&mut ec_aggregate, &protein.ec_numbers);
            Self::count_occurrences(&mut go_aggregate, &protein.go_terms);
            Self::count_occurrences(&mut interpro_aggregate, &protein.interpro);
        }

        Self {
            ec: ec_aggregate,
            go: go_aggregate,
            interpro: interpro_aggregate,
        }
    }

    pub fn create_empty() -> Self {
        Self {
            ec: HashMap::new(),
            go: HashMap::new(),
            interpro: HashMap::new(),
        }
    }
}
