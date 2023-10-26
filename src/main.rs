use clap::Parser;

use rust_implementations::{Arguments, run};

fn main() {
    let args = Arguments::parse();
    run(args);
}
