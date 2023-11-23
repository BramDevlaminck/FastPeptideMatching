use clap::{arg, Parser, ValueEnum};

use tsv_utils::{get_proteins_from_database_file, proteins_to_concatenated_string};

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
    let data = proteins_to_concatenated_string(&proteins);
    let u8_text = data.as_bytes();

    let suffix_array = match args.algorithm {
        Algorithm::DivsufsortC => cdivsufsort::sort(u8_text),
        Algorithm::DivsufsortRust => divsufsort::sort(u8_text)
    };
    let (_text, sa) = suffix_array.into_parts();
    println!("{}", sa.len());
}