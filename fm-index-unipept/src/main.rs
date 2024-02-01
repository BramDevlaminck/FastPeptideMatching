use std::error::Error;
use clap::Parser;
use fm_index_unipept::{Arguments, run};

fn main() -> Result<(), Box<dyn Error>>{
    let args = Arguments::parse();
    run(args)?;
    Ok(())
}
