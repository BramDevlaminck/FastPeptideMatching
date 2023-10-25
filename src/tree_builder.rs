use crate::cursor::{Cursor, CursorIterator};
use crate::tree::{Tree};

pub trait TreeBuilder {
    fn new() -> Self;

    fn build(&self, data: &String, tree: Tree) -> Tree;
}

pub struct NaiveBuilder {}

impl TreeBuilder for NaiveBuilder {
    fn new() -> Self {
        NaiveBuilder {}
    }

    fn build(&self, data: &String, mut tree: Tree) -> Tree {
        let mut cursor = Cursor::new(&mut tree);
        let input_string = data.as_bytes();
        let end_index = input_string.len();
        for (i, character) in data.char_indices() {
            let mut index_in_entry = i;
            let mut ret_value = cursor.next(character, input_string);

            while ret_value == CursorIterator::Ok {
                index_in_entry += 1;
                if index_in_entry == end_index {
                    break;
                }
                ret_value = cursor.next(input_string[index_in_entry] as char, input_string);
            }
            // let current_node = &mut tree.arena[cursor.current_node_index_in_arena];
            if ret_value == CursorIterator::InWord {
                cursor.split_and_add(index_in_entry, end_index, input_string);
            } else {
                cursor.add_leaf(end_index, index_in_entry, input_string);
            }
            cursor.reset();
        }

        tree
    }
}
