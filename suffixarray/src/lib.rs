use std::error::Error;
use std::num::NonZeroUsize;

use bitarray::BitArray;
use clap::{arg, Parser, ValueEnum};

use sa_mappings::functionality::FunctionAggregator;
use sa_mappings::proteins::Proteins;
use sa_mappings::taxonomy::{AggregationMethod, TaxonAggregator};
use suffixarray_builder::{build_sa, SAConstructionAlgorithm};
use suffixarray_builder::binary::{load_binary, write_binary};
use tsv_utils::read_lines;

use crate::peptide_search::search_all_peptides;
use crate::sa_searcher::Searcher;
use crate::suffix_to_protein_index::{
    DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinIndex, SuffixToProteinMappingStyle,
};
use crate::util::get_time_ms;

pub mod peptide_search;
pub mod sa_searcher;
pub mod suffix_to_protein_index;
pub mod util;

/// Enum that represents the 5 kinds of search that we support
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    MinMaxBound,
    AllOccurrences,
    TaxonId,
    Analyses,
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    search_file: Option<String>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
    /// This will only build the tree and stop after that is completed. Used during benchmarking.
    #[arg(long)]
    build_only: bool,
    /// Output file to store the built index.
    #[arg(short, long)]
    output: Option<String>,
    /// The sample rate used on the suffix array (default value 1, which means every value in the SA is used)
    #[arg(long, default_value_t = 1)]
    sample_rate: u8,
    /// Set the style used to map back from the suffix to the protein. 2 options <sparse> or <dense>. Dense is default
    /// Dense uses O(n) memory with n the size of the input text, and takes O(1) time to find the mapping
    /// Sparse uses O(m) memory with m the number of proteins, and takes O(log m) to find the mapping
    #[arg(long, value_enum, default_value_t = SuffixToProteinMappingStyle::Sparse)]
    suffix_to_protein_mapping: SuffixToProteinMappingStyle,
    #[arg(long)]
    load_index: Option<String>,
    #[arg(short, long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    construction_algorithm: SAConstructionAlgorithm,
    /// Assume the resulting taxon ID is root (1) whenever a peptide matches >= cutoff proteins
    #[arg(long, default_value_t = 10000)]
    cutoff: usize,
    #[arg(long)]
    threads: Option<NonZeroUsize>,
    #[arg(long)]
    equalize_i_and_l: bool,
    #[arg(long)]
    clean_taxa: bool,
}

pub fn run(mut args: Arguments) -> Result<(), Box<dyn Error>> {
    let taxon_id_calculator =
        TaxonAggregator::try_from_taxonomy_file(&args.taxonomy, AggregationMethod::LcaStar)?;

    let sa = match &args.load_index {
        // load SA from file
        Some(index_file_name) => {
            let (sample_rate, sa) = load_binary(index_file_name)?;
            args.sample_rate = sample_rate;
            // println!("Loading the SA took {} ms and loading the proteins + SA took {} ms", end_loading_ms - start_loading_ms, end_loading_ms - start_reading_proteins_ms);
            // TODO: some kind of security check that the loaded database file and SA match
            sa
        }
        // build the SA
        None => {
            let protein_sequences =
                Proteins::try_from_database_file(&args.database_file, &taxon_id_calculator)?;
            // build_sa(
            //     &mut protein_sequences.input_string.clone(),
            //     &args.construction_algorithm,
            //     args.sample_rate,
            // )?
            BitArray::<38>::with_capacity(1)
        }
    };

    let proteins = Proteins::try_from_database_file(&args.database_file, &taxon_id_calculator)?;

    // if let Some(output) = &args.output {
    //     write_binary(args.sample_rate, &sa, output)?;
    // }

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return Ok(());
    }

    // build the right mapping index, use box to be able to store both types in this variable
    let suffix_index_to_protein: Box<dyn SuffixToProteinIndex> =
        match args.suffix_to_protein_mapping {
            SuffixToProteinMappingStyle::Dense => {
                Box::new(DenseSuffixToProtein::new(&proteins.input_string))
            }
            SuffixToProteinMappingStyle::Sparse => {
                Box::new(SparseSuffixToProtein::new(&proteins.input_string))
            }
        };

    let functional_aggregator = FunctionAggregator {};

    let searcher = Searcher::new(
        sa,
        args.sample_rate,
        suffix_index_to_protein,
        proteins,
        taxon_id_calculator,
        functional_aggregator,
    );

    execute_search(&searcher, &args)?;
    Ok(())
}

/// Perform the search as set with the commandline arguments
fn execute_search(searcher: &Searcher, args: &Arguments) -> Result<(), Box<dyn Error>> {
    let cutoff = args.cutoff;
    let search_file = args
        .search_file
        .as_ref()
        .ok_or("No peptide file provided to search in the database")?;

    let start_time = get_time_ms()?;
    let lines = read_lines(search_file)?;
    let all_peptides: Vec<String> = lines.map_while(Result::ok).collect();

    // Explicitly set the number of threads to use if the commandline argument was set
    if let Some(threads) = args.threads {
        rayon::ThreadPoolBuilder::new()
            .num_threads(threads.get())
            .build_global()?;
    }

    let search_result = search_all_peptides(
        searcher,
        &all_peptides,
        cutoff,
        args.equalize_i_and_l,
        args.clean_taxa,
    );
    println!("{}", serde_json::to_string(&search_result)?);
    let end_time = get_time_ms()?;

    // output to other channel to prevent integrating it into the actual output
    eprintln!(
        "Spend {} ms to search the whole file",
        end_time - start_time
    );

    Ok(())
}

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    fn is_null(&self) -> bool;
}

impl Nullable<u32> for u32 {
    const NULL: u32 = u32::MAX;

    fn is_null(&self) -> bool {
        *self == Self::NULL
    }
}
