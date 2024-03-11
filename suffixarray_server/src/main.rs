use std::collections::HashMap;
use std::error::Error;
use std::sync::Arc;

use axum::{http::StatusCode, Json, Router};
use axum::extract::{DefaultBodyLimit, State};
use axum::routing::{get, post};
use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use suffixarray::binary::load_binary;
use suffixarray::functional_annotations::{FunctionalAnnotationsCounts, PeptideSearchResult};
use suffixarray::searcher::Searcher;
use suffixarray::suffix_to_protein_index::SparseSuffixToProtein;
use tsv_utils::get_proteins_from_database_file;
use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};

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

#[derive(Debug, Deserialize, Serialize)]
struct InputData {
    peptides: Vec<String>,
    cutoff: Option<usize>
}

#[derive(Debug, Serialize)]
struct OutputData {
    result: Vec<SearchResult>,
}

#[derive(Debug, Serialize)]
struct FunctionalAnnotations {
    counts: FunctionalAnnotationsCounts,
    data: HashMap<String, u64>,
}

#[derive(Debug, Serialize)]
struct SearchResult {
    sequence: String,
    lca: usize,
    fa: FunctionalAnnotations,
    taxa: Vec<usize>,
    uniprot_ids: Vec<String>,
    cutoff_used: bool
}

// basic handler that responds with a static string
async fn root() -> &'static str {
    "Server is online"
}

/// Executes the kind of search indicated by the commandline arguments
pub fn search_peptide(
    searcher: &Searcher,
    word: &str,
    cutoff: usize,
) -> Option<(usize, bool, PeptideSearchResult)> {
    let word = word.to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if word.len() < searcher.sample_rate as usize {
        return None;
    }

    let (cutoff_used, suffixes) = searcher.search_matching_suffixes(word.as_bytes(), cutoff);
    let proteins = searcher.retrieve_proteins(&suffixes);
    let annotations = PeptideSearchResult::new(&proteins);
    if cutoff_used {
        Some((1, cutoff_used, annotations))
    } else {
        let id = searcher.retrieve_taxon_id(&proteins);
        id.map(|id_unwrapped| (id_unwrapped, cutoff_used, annotations))
    }
}

async fn calculate(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData>, StatusCode> {
    
    let cutoff = data.cutoff.unwrap_or(10000);
    
    let res: Vec<SearchResult> = data
        .peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(&searcher, peptide, cutoff))
        .enumerate()
        // remove the peptides that did not match any proteins
        .filter(|(_, search_result)| search_result.is_some())
        // transform search results into output
        .map(|(index, search_results)| {
            let (taxon_id, cutoff_used, pept_results) = search_results.unwrap();
            SearchResult {
                sequence: data.peptides[index].clone(),
                lca: taxon_id,
                fa: FunctionalAnnotations { counts: pept_results.counts, data: pept_results.data },
                taxa: pept_results.taxa,
                uniprot_ids: pept_results.uniprot_ids,
                cutoff_used
            }
        })
        .collect();

    // Create output JSON
    let output_data = OutputData { result: res };

    Ok(Json(output_data))
}

#[tokio::main]
async fn main() -> Result<(), Box<dyn Error>> {
    let args = Arguments::parse();
    let Arguments {
        database_file,
        index_file,
        taxonomy,
    } = args;
    // let (sample_rate, sa) = load_binary(&index_file)?;

    let taxon_id_calculator = TaxonIdCalculator::new(&taxonomy, AggregationMethod::LcaStar);
    let proteins = get_proteins_from_database_file(&database_file, &*taxon_id_calculator)?;
    let suffix_index_to_protein = Box::new(SparseSuffixToProtein::new(&proteins.input_string));

    // let searcher = Arc::new(Searcher::new(
    //     sa,
    //     sample_rate,
    //     suffix_index_to_protein,
    //     proteins,
    //     *taxon_id_calculator,
    // ));
    // 
    // // build our application with a route
    // let app = Router::new()
    //     // `GET /` goes to `root`
    //     .route("/", get(root))
    //     // `POST /search` goes to `calculate` and set max payload size to 5 MB
    //     .route("/search", post(calculate)).layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
    //     .with_state(searcher);
    // 
    // let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;
    // println!("server is ready...");
    // axum::serve(listener, app).await?;

    Ok(())
}
