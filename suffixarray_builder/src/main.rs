use std::error::Error;
use clap::Parser;
use suffixarray_builder::{Arguments, build_sa};
use suffixarray_builder::binary::write_binary;
use tsv_utils::get_text_from_database_file;
use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};

fn main() -> Result<(), Box<dyn Error>>{
    let args = Arguments::parse();
    let Arguments { database_file, taxonomy, output, sample_rate, construction_algorithm } = args;
    let taxon_id_calculator = TaxonIdCalculator::new(&taxonomy, AggregationMethod::LcaStar);
    let data = get_text_from_database_file(&database_file, &*taxon_id_calculator)?;
    let sa = build_sa(&data.into_bytes(), &construction_algorithm, sample_rate)?;
    write_binary(sample_rate, &sa, &output)?;
    Ok(())
}