use std::error::Error;
use std::io;
use std::io::Write;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{arg, Parser, ValueEnum};
use get_size::GetSize;

use tsv_utils::{get_proteins_from_database_file, Proteins, read_lines};
use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};

use crate::binary::{load_binary, write_binary};
use crate::searcher::Searcher;
use crate::suffix_to_protein_index::{DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinIndex, SuffixToProteinMappingStyle};

mod searcher;
mod binary;
mod suffix_to_protein_index;

/// Enum that represents the 2 kinds of search that we support
/// - Search until match and return boolean that indicates if there is a match
/// - Search until match, if there is a match return the min and max index in the SA that matches
/// - Search until match, if there is a match search the whole subtree to find all matching proteins
/// - Search until match, there we can immediately retrieve the taxonId that represents all the children
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    MinMaxBound,
    AllOccurrences,
    TaxonId,
}

#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SAConstructionAlgorithm {
    LibDivSufSort,
    LibSais,
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    search_file: Option<String>,
    /// `match` will only look if there is match.
    /// `all-occurrences` will search for the match and look for all the different matches in the subtree.
    /// `min-max-bound` will search for the match and retrieve the minimum and maximum index in the SA that contains a suffix that matches.
    /// `Taxon-id` will search for the matching taxon id using lca*
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
    /// This will change the output to <found (0 or 1)>;<protein length>;<search time in ms>
    /// The given num will be used to run the search x times and the average of these x runs will be given as search time
    #[arg(short, long)]
    verbose: Option<u8>,
    /// This will only build the tree and stop after that is completed. Used during benchmarking.
    #[arg(long)]
    build_only: bool,
    /// Output file to store the built index.
    #[arg(short, long)]
    output: Option<String>,
    /// The sample rate used on the suffix array (default value 1, which means every value in the SA is used)
    #[arg(long, default_value_t = 1)]
    sample_rate: u8, // TODO: is u8 enough for a sample rate?
    /// Set the style used to map back from the suffix to the protein. 2 options <sparse> or <dense>. Dense is default
    /// Dense uses O(n) memory with n the size of the input text, and takes O(1) time to find the mapping
    /// Sparse uses O(m) memory with m the number of proteins, and takes O(log m) to find the mapping
    #[arg(long, value_enum, default_value_t = SuffixToProteinMappingStyle::Dense)]
    suffix_to_protein_mapping: SuffixToProteinMappingStyle,
    #[arg(long)]
    load_index: Option<String>,
    #[arg(short, long, value_enum, default_value_t = SAConstructionAlgorithm::LibSais)]
    construction_algorithm: SAConstructionAlgorithm
}

pub fn run(mut args: Arguments) -> Result<(), Box<dyn Error>> {
    let start_reading_proteins_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    let taxon_id_calculator = TaxonIdCalculator::new(&args.taxonomy, AggregationMethod::LcaStar);
    // println!("taxonomy calculator built");

    let proteins = get_proteins_from_database_file(&args.database_file, &*taxon_id_calculator);

    // construct the sequence that will be used to build the tree
    // println!("read all proteins");
    let current = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    // println!("Time spent for reading: {}", current - start_reading_proteins_ms);

    let sa = match &args.load_index {
        // load SA from file
        Some(index_file_name) => {
            let start_loading_ms = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards").as_nanos() as f64 * 1e-6;
            let (sample_rate, sa) = load_binary(index_file_name)?;
            args.sample_rate = sample_rate;
            let end_loading_ms = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards").as_nanos() as f64 * 1e-6;
            // println!("Loading the SA took {} ms and loading the proteins + SA took {} ms", end_loading_ms - start_loading_ms, end_loading_ms - start_reading_proteins_ms);
            // TODO: some kind of security check that the loaded database file and SA match
            sa
        },
        // build the SA
        None => {
            let mut sa = match &args.construction_algorithm {
                SAConstructionAlgorithm::LibSais => libsais64_rs::sais64(&proteins.input_string),
                SAConstructionAlgorithm::LibDivSufSort => libdivsufsort_rs::divsufsort64(&proteins.input_string)
            }.ok_or("Building suffix array failed")?;
            // println!("SA constructed");

            // make the SA sparse and decrease the vector size if we have sampling (== sampling_rate > 1)
            if args.sample_rate > 1 {
                let mut current_sampled_index = 0;
                for i in 0..sa.len() {
                    let current_sa_val = sa[i];
                    if current_sa_val % args.sample_rate as i64 == 0 {
                        sa[current_sampled_index] = current_sa_val;
                        current_sampled_index += 1;
                    }
                }
                // make shorter
                sa.resize(current_sampled_index, 0);
                // println!("SA is sparse with sampling factor {}", args.sample_rate);
            }
            sa
        }
    };

    if let Some(output) = &args.output {
        // println!("storing index to file {}", output);
        write_binary(args.sample_rate, &sa, &proteins.input_string, output).unwrap();
        // println!("Index written away");
    }

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return Ok(());
    } else if args.mode.is_none() {
        eprintln!("search mode expected!");
        std::process::exit(1);
    }

    // build the right mapping index, use box to be able to store both types in this variable
    let suffix_index_to_protein: Box<dyn SuffixToProteinIndex> = match args.suffix_to_protein_mapping {
        SuffixToProteinMappingStyle::Dense => Box::new(DenseSuffixToProtein::new(&proteins.input_string)),
        SuffixToProteinMappingStyle::Sparse => Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
    };
    // println!("mapping built");

    let searcher = Searcher::new(&sa, args.sample_rate, suffix_index_to_protein.as_ref(), &proteins, &taxon_id_calculator);
    execute_search(searcher, &proteins, &args);
    Ok(())
}

/// Perform the search as set with the commandline arguments
fn execute_search(mut searcher: Searcher, proteins: &Proteins, args: &Arguments) {
    let mode = args.mode.as_ref().unwrap();
    let verbose = args.verbose;
    let mut verbose_output: Vec<String> = vec![];
    if let Some(search_file) = &args.search_file {
        // File `search_file` must exist in the current path
        if let Ok(lines) = read_lines(search_file) {
            let start_time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards").as_nanos() as f64 * 1e-6;
            for line in lines.into_iter().map_while(Result::ok) {
                handle_search_word(&mut searcher, proteins, line, mode, verbose, &mut verbose_output);
            }
            let end_time = SystemTime::now()
                .duration_since(UNIX_EPOCH)
                .expect("Time went backwards").as_nanos() as f64 * 1e-6;
            eprintln!("Spend {} ms to search the whole file", end_time - start_time);
        } else {
            eprintln!("File {} could not be opened!", search_file);
            std::process::exit(1);
        }
    } else {
        loop {
            print!("Input your search string: ");
            io::stdout().flush().unwrap();
            let mut word = String::new();

            if io::stdin().read_line(&mut word).is_err() {
                continue;
            }
            handle_search_word(&mut searcher, proteins, word, mode, verbose, &mut verbose_output);
        }
    }
    verbose_output.iter().for_each(|val| println!("{}", val));
}

fn time_execution(searcher: &mut Searcher, f: &dyn Fn(&mut Searcher) -> bool) -> (bool, f64) {
    let start_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    let found = f(searcher);
    let end_ms = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
    (found, end_ms - start_ms, )
}


/// Executes the kind of search indicated by the commandline arguments
fn handle_search_word(searcher: &mut Searcher, proteins: &Proteins, word: String, search_mode: &SearchMode, verbose: Option<u8>, verbose_output: &mut Vec<String>) {
    let word = match word.strip_suffix('\n') {
        None => word,
        Some(stripped) => String::from(stripped)
    }.to_uppercase();
    // words that are shorter than the sample rate are not searchable
    if word.len() < searcher.sample_rate as usize {
        if verbose.is_some() {
            verbose_output.push(format!("(word too short short for SA sample size);{};", word.len()));
        } else {
            println!("/ (word too short short for SA sample size)")
        }
        return;
    }

    if let Some(num_iter) = verbose {
        let mut found_total: bool = false;
        let mut total_time: f64 = 0.0;
        for _ in 0..num_iter {
            let (found, execution_time) = match *search_mode {
                SearchMode::Match => time_execution(searcher, &|searcher| searcher.search_if_match(word.as_bytes())),
                SearchMode::MinMaxBound => time_execution(searcher, &|searcher| searcher.search_bounds(word.as_bytes()).0),
                SearchMode::AllOccurrences => time_execution(searcher, &|searcher| !searcher.search_protein(word.as_bytes()).is_empty()),
                SearchMode::TaxonId => time_execution(searcher, &|searcher| searcher.search_taxon_id(word.as_bytes()).is_some()),
            };
            total_time += execution_time;
            found_total = found;
        }
        let avg = total_time / (num_iter as f64);

        verbose_output.push(format!("{};{};{}", found_total as u8, word.len(), avg));
    } else {
        match *search_mode {
            SearchMode::Match => println!("{}", searcher.search_if_match(word.as_bytes())),
            SearchMode::MinMaxBound => {
                let mut found= false;
                let mut min_bound = 0;
                let mut max_bound= 0;
                let mut total_time = 0.0;
                for _ in 0..3 {
                    let start_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    (found, min_bound, max_bound) = searcher.search_bounds(word.as_bytes());
                    let end_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    total_time += end_time - start_time;
                }
                let time_spent_searching = total_time/3.0;
                println!("{found};{min_bound};{max_bound};{time_spent_searching}");
            }
            SearchMode::AllOccurrences => {
                let mut results = vec![];
                let mut total_time = 0.0;
                for _ in 0..3 {
                    let start_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    results = searcher.search_protein(word.as_bytes());

                    let end_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    total_time += end_time - start_time;
                }
                let time_spent_searching = total_time/3.0;
                let number_of_proteins = results.len();
                let peptide_length = word.len();
                println!("{peptide_length};{number_of_proteins};{time_spent_searching}");
            }
            SearchMode::TaxonId => {
                let mut result = None;
                let mut total_time = 0.0;
                for _ in 0..1 {
                    let start_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    let max_matches = 4000;
                    let suffixes = searcher.search_matching_suffixes(word.as_bytes(), max_matches);
                    if suffixes.len() >= max_matches {
                        result = Some(1);
                    } else {
                        let proteins = searcher.retrieve_proteins(&suffixes);
                        result = searcher.retrieve_taxon_id(&proteins);
                    }

                    let end_time = SystemTime::now()
                        .duration_since(UNIX_EPOCH)
                        .expect("Time went backwards").as_nanos() as f64 * 1e-6;
                    total_time += end_time - start_time;
                }
                let taxon_id = if let Some(id) = result {
                    format!("{id}")
                } else {
                    "/".to_string()
                };
                println!("{};{};{}", word.len(), taxon_id, total_time / 1.0);
            }
        }
    }
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
