use std::error::Error;
use std::sync::Arc;

use axum::{http::StatusCode, Json, Router};
use axum::extract::{DefaultBodyLimit, State};
use axum::routing::{get, post};
use clap::Parser;
use serde::{Deserialize, Serialize};

use sa_mappings::functionality::FunctionAggregator;
use sa_mappings::proteins::Proteins;
use sa_mappings::taxonomy::{AggregationMethod, TaxonAggregator};
use suffixarray::peptide_search::{
    OutputData, search_all_peptides,
};
use suffixarray::sa_searcher::Searcher;
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

#[allow(dead_code)]
fn default_true() -> bool {
    true
}

#[derive(Debug, Deserialize, Serialize)]
#[allow(non_snake_case)]
struct InputData {
    peptides: Vec<String>,
    #[serde(default = "default_cutoff")] // default value is 10000
    cutoff: usize,
    #[serde(default = "bool::default")]
    // default value is false // TODO: maybe default should be true?
    equalize_I_and_L: bool,
    #[serde(default = "bool::default")] // default value is false
    clean_taxa: bool,
    #[serde(default = "bool::default")] // default value is false
    search_only: bool,
}

// basic handler that responds with a static string
async fn root() -> &'static str {
    "Server is online"
}

/// Search all the peptides in the given data json using the searcher.
async fn search(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData>, StatusCode> {
    let search_result = search_all_peptides(
        &searcher,
        &data.peptides,
        data.cutoff,
        data.equalize_I_and_L,
        data.clean_taxa,
        data.search_only
    );

    Ok(Json(search_result))
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
    let taxon_id_calculator =
        TaxonAggregator::try_from_taxonomy_file(&taxonomy, AggregationMethod::LcaStar)?;

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
        function_aggregator,
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
