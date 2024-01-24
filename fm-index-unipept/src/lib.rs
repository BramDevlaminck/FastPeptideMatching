use std::error::Error;
use std::fs;
use clap::Parser;
use fm_index::converter::RangeConverter;
use fm_index::FMIndex;
use fm_index::suffix_array::{NullSampler, SuffixOrderSampler};
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

    let u8_text = proteins.input_string.as_bytes().to_vec();

    // Converter converts each character into packed representation.
    // `' '` ~ `'~'` represents a range of ASCII printable characters.
    let converter = RangeConverter::new(b'$', b'Z');

    // To perform locate queries, we need to retain suffix array generated in the construction phase.
    // However, we don't need the whole array since we can interpolate missing elements in a suffix array from others.
    // A sampler will _sieve_ a suffix array for this purpose.
    // You can also use `NullSampler` if you don't perform location queries (disabled in type-level).
    // let sampler =  SuffixOrderSampler::new().level(2);
    let sampler =  NullSampler::new();
    let _index = FMIndex::new(u8_text, converter, sampler);
    
    Ok(())
}

