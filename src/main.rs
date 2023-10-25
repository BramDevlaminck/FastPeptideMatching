mod tree_builder;
mod tree;
mod cursor;
mod searcher;
mod read_only_cursor;

use std::{env, fs, io};
use std::io::Write;
use crate::searcher::Searcher;
use crate::tree::{Node, Tree};
use crate::tree_builder::{UkkonenBuilder, TreeBuilder};


fn main() {
    let args: Vec<String> = env::args().collect();

    let file_name = args.get(1).expect("Input file name expected");

    let mut data = fs::read_to_string(file_name).expect("Reading input file to build tree from went wrong");
    data = data.to_uppercase();
    data.push('$');

    let tree = Tree::new(&data, UkkonenBuilder::new());

    let mut searcher = Searcher::new(&tree, data.as_bytes());

    loop {
        print!("Input your search string: ");
        io::stdout().flush().unwrap();
        let mut word = String::new();

        if io::stdin().read_line(&mut word).is_err() {
            continue;
        }
        word = match word.strip_suffix('\n') {
            None => word,
            Some(stripped) => stripped.to_string()
        }.to_uppercase();
        // println!("{}", searcher.search_if_match(word.as_bytes()))
        let results = searcher.search_protein(word.as_bytes());
        println!("found {} matches", results.len());
        results.iter()
            .for_each(|res| println!("* {}", res));
    }
}
