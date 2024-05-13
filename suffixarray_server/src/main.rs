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
use suffixarray::peptide_search::{OutputData, analyse_all_peptides, SearchResultWithAnalysis, SearchOnlyResult, search_all_peptides};
use suffixarray::sa_searcher::Searcher;
use suffixarray::suffix_to_protein_index::SparseSuffixToProtein;
use suffixarray_builder::binary::load_suffix_array;

/// Enum that represents all possible commandline arguments
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

/// Function used by serde to use `true` as a default value
#[allow(dead_code)]
fn default_true() -> bool {
    true
}

/// Struct representing the input arguments accepted by the endpoints
/// 
/// # Arguments
/// * `peptides` - List of peptides we want to process
/// * `cutoff` - The maximum amount of matches to process, default value 10000
/// * `equalize_I_and_L` - True if we want to equalize I and L during search
/// * `clean_taxa` - True if we only want to use proteins marked as "valid"
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
}

#[tokio::main]
async fn main() {
    let args = Arguments::parse();
    if let Err(err) = start_server(args).await {
        eprintln!("{}", err);
        std::process::exit(1);
    }
}

/// Basic handler used to check the server status
async fn root() -> &'static str {
    "Server is online"
}

/// Endpoint executed for peptide matching and taxonomic and functional analysis
///
/// # Arguments
/// * `state(searcher)` - The searcher object provided by the server
/// * `data` - InputData object provided by the user with the peptides to be searched and the config
/// 
/// # Returns
///
/// Returns the search and analysis results from the index as a JSON
async fn analysis(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData<SearchResultWithAnalysis>>, StatusCode> {
    let search_result = analyse_all_peptides(
        &searcher,
        &data.peptides,
        data.cutoff,
        data.equalize_I_and_L,
        data.clean_taxa,
    );

    Ok(Json(search_result))
}

/// Endpoint executed for peptide matching, without any analysis
///
/// # Arguments
/// * `state(searcher)` - The searcher object provided by the server
/// * `data` - InputData object provided by the user with the peptides to be searched and the config
///
/// # Returns
///
/// Returns the search results from the index as a JSON
async fn search(
    State(searcher): State<Arc<Searcher>>,
    data: Json<InputData>,
) -> Result<Json<OutputData<SearchOnlyResult>>, StatusCode> {
    let search_result = search_all_peptides(
        &searcher,
        &data.peptides,
        data.cutoff,
        data.equalize_I_and_L,
        data.clean_taxa,
    );

    Ok(Json(search_result))
}

/// Starts the server with the provided commandline arguments
///
/// # Arguments
/// * `args` - The provided commandline arguments
///
/// # Returns
///
/// Returns ()
/// 
/// # Errors
/// 
/// Returns any error occurring during the startup or uptime of the server
async fn start_server(args: Arguments) -> Result<(), Box<dyn Error>> {
    let Arguments {
        database_file,
        index_file,
        taxonomy,
    } = args;

    eprintln!("Loading suffix array...");
    let (sparseness_factor, sa) = load_suffix_array(&index_file)?;

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
        sparseness_factor,
        suffix_index_to_protein,
        proteins,
        taxon_id_calculator,
        function_aggregator,
    ));

    // build our application with a route
    let app = Router::new()
        // `GET /` goes to `root`
        .route("/", get(root))
        // `POST /analysis` goes to `analysis` and set max payload size to 5 MB
        .route("/analysis", post(analysis))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher.clone())
        // `POST /search` goes to `search` and set max payload size to 5 MB
        .route("/search", post(search))
        .layer(DefaultBodyLimit::max(5 * 10_usize.pow(6)))
        .with_state(searcher);

    let listener = tokio::net::TcpListener::bind("0.0.0.0:3000").await?;
    println!("server is ready...");
    axum::serve(listener, app).await?;

    Ok(())
}
