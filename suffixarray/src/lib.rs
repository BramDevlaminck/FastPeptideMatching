mod searcher;

use std::io;
use std::io::Write;
use clap::{arg, Parser, ValueEnum};

use tsv_utils::{get_proteins_from_database_file, proteins_to_concatenated_string};
use crate::searcher::Searcher;

#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum Algorithm {
    DivsufsortC,
    DivsufsortRust
}

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
    #[arg(short, long, value_enum)]
    algorithm: Algorithm
}

pub fn run(args: Arguments) {
    let proteins = get_proteins_from_database_file(&args.database_file);
    // construct the sequence that will be used to build the tree
    let mut data = proteins_to_concatenated_string(&proteins);
    data = "BANANA$".to_string(); // TODO: remove en remove mut above this
    let u8_text = data.as_bytes();

    let suffix_array = match args.algorithm {
        Algorithm::DivsufsortC => cdivsufsort::sort(u8_text),
        Algorithm::DivsufsortRust => divsufsort::sort(u8_text)
    };
    let (text, sa) = suffix_array.into_parts();
    println!("{}", sa.len());
    println!("{}", sa[0]);

    let mut current_protein_index: u32 = 0;
    let mut index_to_protein: Vec<Option<u32>> = vec![];
    for &char in text.iter() {
        if char == b'%' || char == b'$' {
            current_protein_index += 1;
            index_to_protein.push(None);
        } else {
            index_to_protein.push(Some(current_protein_index));
        }
    }

    let searcher = Searcher::new(text, &sa, &proteins);
    loop {
        print!("Input your search string: ");
        io::stdout().flush().unwrap();
        let mut word = String::new();

        if io::stdin().read_line(&mut word).is_err() {
            continue;
        }
        word = match word.strip_suffix('\n') {
            None => word,
            Some(stripped) => String::from(stripped)
        }.to_uppercase();
        let out = searcher.search_bounds(word.as_bytes());
        println!("{:?}", out);
    }
    // search

}