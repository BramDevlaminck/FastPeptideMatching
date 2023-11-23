use clap::Parser;
use suffixarray::{Arguments, run};

fn main() {
    let args = Arguments::parse();
    run(args);
}
