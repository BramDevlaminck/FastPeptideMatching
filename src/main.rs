mod tree_builder;
mod tree;
mod cursor;

use std::{env, fs, io};
use std::io::Write;
use crate::tree::{Node, Tree};
use crate::tree_builder::{NaiveBuilder, TreeBuilder};


fn main() {
    let args: Vec<String> = env::args().collect();

    let file_name = args.get(1).expect("Input file name expected");

    let mut data = fs::read_to_string(file_name).expect("Reading input file to build tree from went wrong");
    data = data.to_uppercase();
    data.push('$');

    let tree = Tree::new(data, NaiveBuilder::new());

    println!("{:?}", tree);

    loop {
        print!("Input your search string: ");
        io::stdout().flush().unwrap();
        let mut word = String::new();

        if io::stdin().read_line(&mut word).is_err() {
            continue;
        }

        println!("{}", word);
    }
}
