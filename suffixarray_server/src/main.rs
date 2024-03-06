use std::error::Error;
use std::sync::Arc;

use axum::{http::StatusCode, Json, Router};
use axum::extract::{DefaultBodyLimit, State};
use axum::routing::{get, post};
use clap::Parser;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use suffixarray::binary::load_binary;
use suffixarray::functional_annotations::FunctionalAnnotations;
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
    // Define your input JSON structure
    peptides: Vec<String>,
}

#[derive(Debug, Serialize)]
struct OutputData {
    // Define your output JSON structure
    result: Vec<SearchResult>,
}

#[derive(Debug, Serialize)]
struct SearchResult {
    peptide: String,
    taxon_id: usize,
    proteins: Vec<String>,
    functional_annotations: FunctionalAnnotations,
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
) -> Option<(usize, Vec<String>, FunctionalAnnotations)> {
    let word = word.to_uppercase();

    // words that are shorter than the sample rate are not searchable
    if word.len() < searcher.sample_rate as usize {
        return None;
    }

    let suffixes = searcher.search_matching_suffixes(word.as_bytes(), cutoff);

    if suffixes.len() >= cutoff {
        Some((1, Vec::new(), FunctionalAnnotations::create_empty()))
    } else {
        let proteins = searcher.retrieve_proteins(&suffixes);
        let id = searcher.retrieve_taxon_id(&proteins);
        if let Some(id_unwrapped) = id {
            let uniprot_ids = FunctionalAnnotations::get_uniprot_ids(&proteins);
            let annotations = FunctionalAnnotations::new(&proteins);
            Some((id_unwrapped, uniprot_ids, annotations))
        } else {
            None
        }
    }
}

async fn calculate(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData>, StatusCode> {
    let res: Vec<SearchResult> = data
        .peptides
        .par_iter()
        // calculate the results
        .map(|peptide| search_peptide(&searcher, peptide, 10000))
        .enumerate()
        // remove the peptides that did not match any proteins
        .filter(|(_, search_result)| search_result.is_some())
        // transform search results into output
        .map(|(index, search_results)| {
            let (taxon_id, proteins, functional_annotations) = search_results.unwrap();
            SearchResult {
                peptide: data.peptides[index].clone(),
                taxon_id,
                proteins,
                functional_annotations,
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
    let (sample_rate, sa) = load_binary(&index_file)?;

    let taxon_id_calculator = TaxonIdCalculator::new(&taxonomy, AggregationMethod::LcaStar);
    let proteins = get_proteins_from_database_file(&database_file, &*taxon_id_calculator)?;
    let suffix_index_to_protein = Box::new(SparseSuffixToProtein::new(&proteins.input_string));

    let searcher = Arc::new(Searcher::new(
        sa,
        sample_rate,
        suffix_index_to_protein,
        proteins,
        *taxon_id_calculator,
    ));
    
    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `POST /search` goes to `calculate` and set max payload size to 5 MB
        .route("/search", post(calculate)).layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;
    axum::serve(listener, app).await?;
    println!("server is ready...");

    Ok(())
}
