use std::error::Error;
use std::sync::Arc;

use axum::extract::{DefaultBodyLimit, State};
use axum::routing::{get, post};
use axum::{http::StatusCode, Json, Router};
use clap::Parser;
use rayon::prelude::*;
use sa_mappings::functionality::{FunctionAggregator, FunctionalAggregation};
use sa_mappings::proteins::Proteins;
use sa_mappings::taxonomy::{AggregationMethod, TaxonAggregator};
use serde::{Deserialize, Serialize};

use suffixarray::searcher::{SearchAllSuffixesResult, Searcher};
use suffixarray::suffix_to_protein_index::SparseSuffixToProtein;
use suffixarray_builder::binary::load_binary;

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    index_file: String,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
}

/// Function used by serde to place a default value in the cutoff field of the input
fn default_cutoff() -> usize {
    10000
}

#[derive(Debug, Deserialize, Serialize)]
#[allow(non_snake_case)]
struct InputData {
    peptides: Vec<String>,
    #[serde(default = "default_cutoff")] // default value is 10000
    cutoff: usize,
    #[serde(default = "bool::default")] // default value is false % TODO: maybe default should be true?
    equalize_I_and_L: bool,
}

#[derive(Debug, Serialize)]
struct OutputData {
    result: Vec<SearchResult>,
}

#[derive(Debug, Serialize)]
struct SearchResult {
    sequence: String,
    lca: Option<usize>,
    taxa: Vec<usize>,
    uniprot_accessions: Vec<String>,
    fa: Option<FunctionalAggregation>,
    cutoff_used: bool,
}

// basic handler that responds with a static string
async fn root() -> &'static str {
    "Server is online"
}

pub struct PeptideSearchResult {
    lca: Option<usize>,
    cutoff_used: bool,
    uniprot_accession_numbers: Vec<String>,
    taxa: Vec<usize>,
    fa: Option<FunctionalAggregation>,
}

/// Executes the search of 1 peptide
pub fn search_peptide(
    searcher: &Searcher,
    peptide: &str,
    cutoff: usize,
    equalize_i_and_l: bool,
) -> Option<PeptideSearchResult> {
    let peptide = peptide.to_uppercase();

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
        SearchAllSuffixesResult::SearchResult(matched_suffixes) => {
            matched_suffixes
        }
        SearchAllSuffixesResult::NoMatches => {
            return None;
        },
    };

    let proteins = searcher.retrieve_proteins(&suffixes);
    let (uniprot_accession_numbers, taxa) = Searcher::get_uniprot_and_taxa_ids(&proteins);
    // calculate the lca
    let lca = if cutoff_used {
        Some(1)
    } else {
        searcher.retrieve_lca(&proteins)
    };
    let fa = searcher.retrieve_function(&proteins);
    // output the result
    Some(PeptideSearchResult {
        lca,
        cutoff_used,
        uniprot_accession_numbers,
        taxa,
        fa
    })
}

/// Search all the peptides in the given data json using the searcher.
async fn search(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData>, StatusCode> {
    let res: Vec<SearchResult> = data
        .peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(&searcher, peptide, data.cutoff, data.equalize_I_and_L))
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
                sequence: data.peptides[index].clone(),
                lca,
                taxa,
                uniprot_accessions: uniprot_accession_numbers,
                fa,
                cutoff_used,
            }
        })
        .collect();

    // Create output JSON
    let output_data = OutputData { result: res };

    Ok(Json(output_data))
}

#[tokio::main]
async fn main() {
    let args = Arguments::parse();
    if let Err(err) = start_server(args).await {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}

async fn start_server(args: Arguments) -> Result<(), Box<dyn Error>> {
    let Arguments {
        database_file,
        index_file,
        taxonomy,
    } = args;

    eprintln!("Loading suffix array...");
    let (sample_rate, sa) = load_binary(&index_file)?;

    eprintln!("Loading taxon file...");
    let taxon_id_calculator = TaxonAggregator::try_from_taxonomy_file(&taxonomy, AggregationMethod::LcaStar)?;

    let function_aggregator = FunctionAggregator {};

    eprintln!("Loading proteins...");
    let proteins = Proteins::try_from_database_file(&database_file, &taxon_id_calculator)?;
    let suffix_index_to_protein = Box::new(SparseSuffixToProtein::new(&proteins.input_string));

    eprintln!("Creating searcher...");
    let searcher = Arc::new(Searcher::new(
        sa,
        sample_rate,
        suffix_index_to_protein,
        proteins,
        taxon_id_calculator,
        function_aggregator
    ));

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `POST /search` goes to `calculate` and set max payload size to 5 MB
        .route("/search", post(search))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;
    println!("server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}
