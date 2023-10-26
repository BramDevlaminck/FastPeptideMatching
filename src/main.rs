mod tree_builder;
mod tree;
mod cursor;
mod searcher;
mod read_only_cursor;
use clap::{Parser, ValueEnum};
use std::{fs, io};
use std::fs::File;
use std::io::{BufRead, Write};
use std::path::Path;
use crate::searcher::Searcher;
use crate::tree::{Node, Tree};
use crate::tree_builder::{UkkonenBuilder, TreeBuilder};
use std::time::{SystemTime, UNIX_EPOCH};

/// Enum that represents the 2 kinds of search that we support
/// - Search until match and return boolean that indicates if there is a match
/// - Search until match, if there is a match search the whole subtree to find all matching proteins
#[derive(ValueEnum, Clone, Debug, PartialEq)]
enum SearchMode {
    Match,
    AllOccurrences
}

#[derive(Parser, Debug)]
struct Arguments {
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long)]
    search_file: Option<String>,
    #[arg(long)]
    build_only: bool,
    #[arg(short, long, value_enum)]
    mode: Option<SearchMode>,
    // this will change the output to <found>;<protein length>;<search time in ms>
    #[arg(short, long, value_enum)]
    verbose: bool
}

// The output is wrapped in a Result to allow matching on errors
// Returns an Iterator to the Reader of the lines of the file.
fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}

fn time_execution(searcher: &mut Searcher,f: &dyn Fn(&mut Searcher) -> bool) -> (bool, f64) {
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
fn handle_search_word(searcher: &mut Searcher, word: String, search_mode: &SearchMode, verbose: bool, verbose_output: &mut Vec<String>) {
    let word = match word.strip_suffix('\n') {
        None => word,
        Some(stripped) => String::from(stripped)
    }.to_uppercase();
    if verbose {
        // initialization with default value needed
        let mut execution_time: f64 = 0.0;
        let mut found: bool = false;

        if *search_mode == SearchMode::Match {
            (found, execution_time) = time_execution(searcher ,& |searcher| searcher.search_if_match(word.as_bytes()));
        } else {
            (found, execution_time) = time_execution(searcher, &|searcher| !searcher.find_all_suffix_indices(word.as_bytes()).is_empty());
        }

        verbose_output.push(format!("{};{};{}", found as u8, word.len(), execution_time));
    } else if *search_mode == SearchMode::Match {
        println!("{}", searcher.search_if_match(word.as_bytes()))
    } else {
        let results = searcher.search_protein(word.as_bytes());
        println!("found {} matches", results.len());
        results.iter()
            .for_each(|res| println!("* {}", res));
    }

}


fn main() {
    let args = Arguments::parse();
    let mut data = fs::read_to_string(args.database_file).expect("Reading input file to build tree from went wrong");
    data = data.to_uppercase();
    data.push('$');

    let tree = Tree::new(&data, UkkonenBuilder::new());

    // option that only builds the tree, but does not allow for querying (easy for benchmark purposes)
    if args.build_only {
        return
    } else if args.mode.is_none() {
        eprintln!("search mode expected!");
        std::process::exit(1);
    }

    let mut searcher = Searcher::new(&tree, data.as_bytes());
    let mode = &args.mode.unwrap();
    let verbose = args.verbose;
    let mut verbose_output: Vec<String> = vec![];
    if let Some(search_file) = args.search_file {
        // File `search_file` must exist in the current path
        if let Ok(lines) = read_lines(&search_file) {
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
