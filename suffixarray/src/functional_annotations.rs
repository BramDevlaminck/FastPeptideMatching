use serde::Serialize;
use std::collections::HashMap;
use umgap::taxon::TaxonId;
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

#[derive(Debug, Default)]
pub struct PeptideSearchResult {
    pub counts: FunctionalAnnotationsCounts,
    pub data: HashMap<String, u64>,
    pub uniprot_ids: Vec<String>,
    pub taxa: Vec<TaxonId>
}

impl PeptideSearchResult {

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
        let mut taxa = vec![];

        for &protein in proteins {
            Self::count_occurrences(&mut annotations_aggregate, &protein.annotations); // TODO: possibly could precalculate this per protein, but will need extra memory
            // count the number of proteins that have an annotation
            let mut has_annotation = false;
            for annotation in &protein.annotations {
                if annotation.starts_with("GO") {
                    go_counts += 1;
                    has_annotation = true;
                } else if annotation.starts_with("IPR") {
                    ipr_counts += 1;
                    has_annotation = true;
                } else { // TODO: when tsv input is fixed, we should do `else if annotation.starts_with("EC")`
                    ec_counts += 1;
                    has_annotation = true;
                }
            }
            
            if has_annotation {
                all_counts += 1;
            }
            
            taxa.push(protein.id);
            uniprot_ids.push(protein.uniprot_id.clone());
        }

        Self {
            counts: FunctionalAnnotationsCounts::new(all_counts, ec_counts, go_counts, ipr_counts),
            data: annotations_aggregate,
            uniprot_ids,
            taxa
        }
    }
}
