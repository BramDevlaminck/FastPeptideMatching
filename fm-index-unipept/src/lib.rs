use std::error::Error;
use clap::Parser;
use tsv_utils::get_proteins_from_database_file;

#[derive(Parser, Debug)]
pub struct Arguments {
    /// File with the proteins used to build the suffix tree. All the proteins are expected to be concatenated using a `#`.
    #[arg(short, long)]
    database_file: String,
}

pub fn run(args: Arguments) -> Result<(), Box<dyn Error>> {
    let proteins = get_proteins_from_database_file(&args.database_file);
    // construct the sequence that will be used to build the tree

    let u8_text = proteins.input_string.as_bytes();
    
    Ok(())
}

