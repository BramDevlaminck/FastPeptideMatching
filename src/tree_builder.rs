use crate::cursor::{Cursor, CursorIterator};
use crate::tree::{Node, Range, Tree};

pub trait TreeBuilder {
    fn new() -> Self;

    fn build(&self, tree: Tree) -> Tree;
}

pub struct NaiveBuilder {}

impl TreeBuilder for NaiveBuilder {
    fn new() -> Self {
        NaiveBuilder {}
    }

    fn build(&self, mut tree: Tree) -> Tree {
        let mut cursor = Cursor::new(&mut tree);
        let input_string = cursor.tree.string.as_bytes();
        let end_index = input_string.len();
        for (i, character) in cursor.tree.string.char_indices() {
            let mut index_in_entry = i;
            let mut ret_value = cursor.next(character);

            while ret_value == CursorIterator::Ok {
                index_in_entry += 1;
                if index_in_entry == end_index {
                    break;
                }
                ret_value = cursor.next(input_string[index_in_entry] as char);
            }

            if ret_value == CursorIterator::InWord {
                let (new_node_index, new_node) = Node::build_in_arena(Range::new(index_in_entry, end_index), cursor.current_node_index_in_arena, std::array::from_fn(|_| None), None, None, &mut cursor.tree);
                let (node_to_insert_in_edge_index, node_to_insert_in_edge) = Node::build_in_arena(Range::new(cursor.current_node.range.start + cursor.index, cursor.current_node.range.end), cursor.current_node_index_in_arena, cursor.current_node.children, None, None, &mut cursor.tree);
                cursor.current_node.set_new_children(
                    vec![
                        (input_string[node_to_insert_in_edge.range.start] as char, node_to_insert_in_edge_index),
                        (input_string[new_node.range.start] as char, new_node_index)
                    ]
                );
                cursor.current_node.range.end = cursor.current_node.range.start + cursor.index;
            } else {
                let (new_node_index, new_node) = Node::build_in_arena(Range::new(index_in_entry, end_index), cursor.current_node_index_in_arena, std::array::from_fn(|_| None), None, None, &mut cursor.tree);
                cursor.current_node.add_child(input_string[index_in_entry] as char, new_node_index);
            }

        }

        tree
    }
}
