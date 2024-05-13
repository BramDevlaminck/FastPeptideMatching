use crate::sa_searcher::{SearchAllSuffixesResult, Searcher};
use rayon::prelude::*;
use sa_mappings::functionality::FunctionalAggregation;
use sa_mappings::proteins::Protein;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct OutputData<T: Serialize> {
    result: Vec<T>,
}

#[derive(Debug, Serialize)]
pub struct SearchResultWithAnalysis {
    sequence: String,
    lca: Option<usize>,
    taxa: Vec<usize>,
    uniprot_accession_numbers: Vec<String>,
    fa: Option<FunctionalAggregation>,
    cutoff_used: bool,
}

#[derive(Debug, Serialize)]
pub struct SearchOnlyResult {
    sequence: String,
    proteins: Vec<ProteinInfo>,
    cutoff_used: bool,
}

#[derive(Debug, Serialize)]
pub struct ProteinInfo {
    taxon: usize,
    uniprot_accession: String,
    functional_annotations: Vec<String>,
}

pub fn analyse_all_peptides(
    searcher: &Searcher,
    peptides: &Vec<String>,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> OutputData<SearchResultWithAnalysis> {
    let res: Vec<SearchResultWithAnalysis> = peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_analysis(searcher, peptide, cutoff, equalize_i_and_l, clean_taxa))
        // remove the peptides that did not match any proteins
        .filter_map(|search_result| search_result)
        .collect();

    OutputData { result: res }
}

pub fn search_all_peptides(
    searcher: &Searcher,
    peptides: &Vec<String>,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> OutputData<SearchOnlyResult> {
    let res: Vec<SearchOnlyResult> = peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide_only(searcher, peptide, cutoff, equalize_i_and_l, clean_taxa))
        // remove the peptides that did not match any proteins
        .filter_map(|search_result| search_result)
        .collect();

    OutputData { result: res }
}

/// Executes the search of 1 peptide
pub fn search_analysis(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> Option<SearchResultWithAnalysis> {
    let (cutoff_used, mut proteins) =
        search_peptide(searcher, peptide, cutoff, equalize_i_and_l, clean_taxa)?;

    if clean_taxa {
        proteins.retain(|protein| searcher.taxon_valid(protein))
    }

    // calculate the lca
    let lca = if cutoff_used {
        Some(1)
    } else {
        searcher.retrieve_lca(&proteins)
    };

    // return None if the LCA is none
    lca?;

    let mut uniprot_accession_numbers = vec![];
    let mut taxa = vec![];

    for protein in &proteins {
        taxa.push(protein.taxon_id);
        uniprot_accession_numbers.push(protein.uniprot_id.clone());
    }

    let fa = searcher.retrieve_function(&proteins);
    // output the result
    Some(SearchResultWithAnalysis {
        sequence: peptide.to_string(),
        lca,
        cutoff_used,
        uniprot_accession_numbers,
        taxa,
        fa,
    })
}

/// Executes the search of 1 peptide
pub fn search_peptide<'a>(
    searcher: &'a Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> Option<(bool, Vec<&'a Protein>)> {
    let peptide = peptide.strip_suffix('\n').unwrap_or(peptide).to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if peptide.len() < searcher.sample_rate as usize {
        return None;
    }

    let suffix_search =
        searcher.search_matching_suffixes(peptide.as_bytes(), cutoff, equalize_i_and_l);
    let mut cutoff_used = false;
    let suffixes = match suffix_search {
        SearchAllSuffixesResult::MaxMatches(matched_suffixes) => {
            cutoff_used = true;
            matched_suffixes
        }
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => matched_suffixes,
        SearchAllSuffixesResult::NoMatches => {
            return None;
        }
    };

    let mut proteins = searcher.retrieve_proteins(&suffixes);
    if clean_taxa {
        proteins.retain(|protein| searcher.taxon_valid(protein))
    }

    Some((cutoff_used, proteins))
}

/// Executes the search of 1 peptide
pub fn search_peptide_only(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> Option<SearchOnlyResult> {
    let (cutoff_used, proteins) =
        search_peptide(searcher, peptide, cutoff, equalize_i_and_l, clean_taxa)?;

    let annotations = searcher.get_all_functional_annotations(&proteins);

    let mut protein_info: Vec<ProteinInfo> = vec![];
    for (&protein, annotations) in proteins.iter().zip(annotations) {
        protein_info.push(ProteinInfo {
            taxon: protein.taxon_id,
            uniprot_accession: protein.uniprot_id.clone(),
            functional_annotations: annotations,
        })
    }

    Some(SearchOnlyResult {
        sequence: peptide.to_string(),
        proteins: protein_info,
        cutoff_used,
    })
}
