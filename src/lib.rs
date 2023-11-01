use std::fs::File;
use std::io;
use std::io::{BufRead, Write};
use std::ops::Add;
use std::path::Path;
use std::time::{SystemTime, UNIX_EPOCH};

use clap::{arg, Parser, ValueEnum};
use umgap::taxon::TaxonId;
use crate::taxon_id_calculator::{TaxonIdCalculator};

use crate::searcher::Searcher;
use crate::tree::Tree;
use crate::tree_builder::{TreeBuilder, UkkonenBuilder};

mod tree_builder;
mod tree;
mod cursor;
mod read_only_cursor;
mod searcher;
mod taxon_id_calculator;


/// Enum that represents the 2 kinds of search that we support
/// - Search until match and return boolean that indicates if there is a match
/// - Search until match, if there is a match search the whole subtree to find all matching proteins
/// - Search until match, there we can immediately retrieve the taxonId that represents all the children
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SearchMode {
    Match,
    AllOccurrences,
    TaxonId,
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    /// A file that contains sequences that we want to search in the tree. Every line contains a new sequence.
    #[arg(short, long)]
    search_file: Option<String>,
    /// This will only build the tree and stop after that is completed. Used during benchmarking.
    #[arg(long)]
    build_only: bool,
    /// `match` will only look if there is match.
    /// While `all-occurrences` will search for the match and look for all the different matches in the subtree.
    /// `Taxon-id` will search for the matching taxon id using lca*
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    /// This will change the output to <found (0 or 1)>;<protein length>;<search time in ms>
    /// The given num will be used to run the search x times and the average of these x runs will be given as search time
    #[arg(short, long)]
    verbose: Option<u8>,
    #[arg(short, long)]
    /// The taxonomy to be used as a tsv file. This is a preprocessed version of the NCBI taxonomy.
    taxonomy: String,
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
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
fn handle_search_word(searcher: &mut Searcher, word: String, search_mode: &SearchMode, verbose: Option<u8>, verbose_output: &mut Vec<String>) {
    let word = match word.strip_suffix('\n') {
        None => word,
        Some(stripped) => String::from(stripped)
    }.to_uppercase();
    if let Some(num_iter) = verbose {
        let mut found_total: bool = false;
        let mut total_time: f64 = 0.0;
        for _ in 0..num_iter {
            let execution_time: f64;
            let found: bool;

            // CLion / RustRover complains about the destructuring into existing variables not working, but this does indeed work (since rust 1.59)
            match *search_mode {
                SearchMode::Match => (found, execution_time) = time_execution(searcher, &|searcher| searcher.search_if_match(word.as_bytes())),
                SearchMode::AllOccurrences => (found, execution_time) = time_execution(searcher, &|searcher| !searcher.find_all_suffix_indices(word.as_bytes()).is_empty()),
                SearchMode::TaxonId => (found, execution_time) = time_execution(searcher, &|searcher| searcher.search_taxon_id(word.as_bytes()).is_some()),
            }
            total_time += execution_time;
            found_total = found;
        }
        let avg = total_time / (num_iter as f64);

        verbose_output.push(format!("{};{};{}", found_total as u8, word.len(), avg));
    } else {
        match *search_mode {
            SearchMode::Match => println!("{}", searcher.search_if_match(word.as_bytes())),
            SearchMode::AllOccurrences => {
                let results = searcher.search_protein(word.as_bytes());
                println!("found {} matches", results.len());
                results.iter()
                    .for_each(|res| println!("* {}", res.sequence));
            }
            SearchMode::TaxonId => {
                match searcher.search_taxon_id(word.as_bytes()) {
                    Some(taxon_id) => println!("{}", taxon_id),
                    None => println!("/"),
                }
            }
        }
    }
}

pub struct Protein {
    pub sequence: String,
    pub id: TaxonId,
}

/// Main run function that executes all the logic with the received arguments
pub fn run(args: Arguments) {
    let mut proteins: Vec<Protein> = vec![];
    if let Ok(lines) = read_lines(&args.database_file) {
        for line in lines.into_iter().flatten() {
            let [_, _, protein_id_str, _, _, protein_sequence]: [&str; 6] = line.splitn(6, '\t').collect::<Vec<&str>>().try_into().unwrap();
            let protein_id_usize = protein_id_str.parse::<TaxonId>().expect("Could not parse id of protein to usize!");
            proteins.push(
                Protein {
                    sequence: protein_sequence.to_uppercase(),
                    id: protein_id_usize,
                }
            )
        }
    } else {
        eprintln!("Database file {} could not be opened!", args.database_file);
        std::process::exit(1);
    }
    let data = proteins
        .iter()
        .map(|prot| prot.sequence.clone())
        .collect::<Vec<String>>()
        .join("#")
        .add("$");

    // build the tree
    let mut tree = Tree::new(&data, UkkonenBuilder::new());
    // fill in the Taxon Ids in the tree using the LCA implementations from UMGAP
    TaxonIdCalculator::new(&args.taxonomy).calculate_taxon_ids_recursive(&mut tree, &proteins);

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return;
    } else if args.mode.is_none() {
        eprintln!("search mode expected!");
        std::process::exit(1);
    }

    let mut searcher = Searcher::new(&tree, data.as_bytes(), &proteins);
    let mode = &args.mode.unwrap();
    let verbose = args.verbose;
    let mut verbose_output: Vec<String> = vec![];
    if let Some(search_file) = &args.search_file {
        // File `search_file` must exist in the current path
        if let Ok(lines) = read_lines(search_file) {
            for line in lines.into_iter().flatten() {
                handle_search_word(&mut searcher, line, mode, verbose, &mut verbose_output);
            }
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
            handle_search_word(&mut searcher, word, mode, verbose, &mut verbose_output);
        }
    }
    verbose_output.iter().for_each(|val| println!("{}", val));
}