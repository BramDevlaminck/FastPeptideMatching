use serde::Serialize;
use std::collections::HashMap;
use tsv_utils::Protein;

#[derive(Debug, Serialize, Default)]
#[allow(non_snake_case)]
pub struct FunctionalAnnotationsCounts {
    all: u64,
    EC: u64,
    GO: u64,
    IPR: u64
}

impl FunctionalAnnotationsCounts {
    pub fn new(all_counts: u64, ec_counts: u64, go_counts: u64, ipr_counts: u64) -> Self {
        Self {
            all: all_counts,
            EC: ec_counts,
            GO: go_counts,
            IPR: ipr_counts
        }
    }
} 

#[derive(Debug, Serialize, Default)]
pub struct FunctionalAnnotations {
    counts: FunctionalAnnotationsCounts,
    data: HashMap<String, u64>,
    uniprot_ids: Vec<String>
}

impl FunctionalAnnotations {

    fn count_occurrences(mapping: &mut HashMap<String, u64>, values: &Vec<String>) {
        for val in values {
            mapping
                .entry(val.clone())
                .and_modify(|v| *v += 1)
                .or_insert(1);
        }
    }

    pub fn new(proteins: &[&Protein]) -> Self {
        let mut annotations_aggregate = HashMap::new();
        let mut all_counts: u64 = 0;
        let mut ec_counts: u64 = 0;
        let mut go_counts: u64 = 0;
        let mut ipr_counts: u64 = 0;

        let mut uniprot_ids = vec![];

        for &protein in proteins {
            Self::count_occurrences(&mut annotations_aggregate, &protein.ec_numbers);
            Self::count_occurrences(&mut annotations_aggregate, &protein.go_terms);
            Self::count_occurrences(&mut annotations_aggregate, &protein.interpro);
            
            // count the number of proteins that have an annotation
            let mut has_annotation = false;
            if !protein.ec_numbers.is_empty() {
                ec_counts += 1;
                has_annotation = true;
            }
            if !protein.go_terms.is_empty() {
                go_counts += 1;
                has_annotation = true;
            }
            if !protein.interpro.is_empty() {
                ipr_counts += 1;
                has_annotation = true;
            }
            if has_annotation {
                all_counts += 1;
            }
            
            uniprot_ids.push(protein.uniprot_id.clone());
        }
        
        Self {
            counts: FunctionalAnnotationsCounts::new(all_counts, ec_counts, go_counts, ipr_counts),
            data: annotations_aggregate,
            uniprot_ids
        }
    }
}
