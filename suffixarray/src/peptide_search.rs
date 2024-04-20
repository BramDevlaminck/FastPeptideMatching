use crate::sa_searcher::{SearchAllSuffixesResult, Searcher};
use rayon::prelude::*;
use sa_mappings::functionality::FunctionalAggregation;
use serde::Serialize;

#[derive(Debug, Serialize)]
pub struct OutputData {
    result: Vec<SearchResult>,
}

#[derive(Debug, Serialize)]
pub struct SearchResult {
    sequence: String,
    lca: Option<usize>,
    taxa: Vec<usize>,
    uniprot_accessions: Vec<String>,
    fa: Option<FunctionalAggregation>,
    cutoff_used: bool,
}

pub struct PeptideSearchResult {
    lca: Option<usize>,
    cutoff_used: bool,
    uniprot_accession_numbers: Vec<String>,
    taxa: Vec<usize>,
    fa: Option<FunctionalAggregation>,
}

pub fn search_all_peptides(
    searcher: &Searcher,
    peptides: &Vec<String>,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> OutputData {
    let res: Vec<SearchResult> = peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(searcher, peptide, cutoff, equalize_i_and_l, clean_taxa))
        .enumerate()
        // remove the peptides that did not match any proteins
        .filter(|(_, search_result)| search_result.is_some())
        // transform search results into output
        .map(|(index, search_results)| {
            let PeptideSearchResult {
                lca,
                cutoff_used,
                uniprot_accession_numbers,
                taxa,
                fa,
            } = search_results.unwrap();
            SearchResult {
                sequence: peptides[index].clone(),
                lca,
                taxa,
                uniprot_accessions: uniprot_accession_numbers,
                fa,
                cutoff_used,
            }
        })
        .collect();

    OutputData { result: res }
}

/// Executes the search of 1 peptide
pub fn search_peptide(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool,
    clean_taxa: bool,
) -> Option<PeptideSearchResult> {
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
    
    let (uniprot_accession_numbers, taxa) = Searcher::get_uniprot_and_taxa_ids(&proteins);
    // calculate the lca
    let lca = if cutoff_used {
        Some(1)
    } else {
        searcher.retrieve_lca(&proteins)
    };
    
    // return None if the LCA is none
    lca?;
    
    let fa = searcher.retrieve_function(&proteins);
    // output the result
    Some(PeptideSearchResult {
        lca,
        cutoff_used,
        uniprot_accession_numbers,
        taxa,
        fa,
    })
}
