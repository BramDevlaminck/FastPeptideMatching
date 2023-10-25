use crate::Node;
use crate::tree::{Range, Tree};

#[derive(Debug, PartialEq)]
pub enum CursorIterator {
    Ok,
    AtEnd,
    InWord
}

pub struct Cursor<'a> {
    pub current_node_index_in_arena: usize,
    // pub current_node: &'a mut Node,
    pub index: usize,
    pub tree: &'a mut Tree
}

impl<'a> Cursor<'a> {

    pub fn new(tree: &'a mut Tree) -> Cursor<'a> {
        Cursor {
            current_node_index_in_arena: 0,
            // current_node: &mut tree.arena[0],
            index: 0,
            tree
        }
    }

    pub fn next(&mut self, next_character: char, bytes_input: &[u8]) -> CursorIterator {
        let current_node = &self.tree.arena[self.current_node_index_in_arena];
        if self.index < current_node.range.length() {
            if bytes_input[current_node.range.start + self.index] as char == next_character {
                self.index += 1;
                return CursorIterator::Ok;
            }
            return CursorIterator::InWord
        }

        if let Some(child) = current_node.get_child(next_character) {
            self.current_node_index_in_arena = child;
            // self.current_node = &mut self.tree.arena[self.current_node_index_in_arena];
            self.index = 1;
            return CursorIterator::Ok
        }

        return CursorIterator::AtEnd
    }

    pub fn reset(&mut self) {
        self.index = 0;
        self.current_node_index_in_arena = 0;
        // self.current_node = &mut self.tree.arena[0];
    }

    pub fn split_and_add(&mut self, index_in_entry: usize, end_index: usize, input_string: &[u8]) {
        let new_node = Node::build_in_arena(Range::new(index_in_entry, end_index), self.current_node_index_in_arena, std::array::from_fn(|_| None), None, None);
        let new_node_char = input_string[new_node.range.start] as char;
        let new_node_index = self.tree.arena.len();
        self.tree.arena.push(new_node);

        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        let node_to_insert_in_edge = Node::build_in_arena(Range::new(current_node.range.start + self.index, current_node.range.end), self.current_node_index_in_arena, current_node.children, None, None);
        current_node.set_new_children(
            vec![
                (input_string[node_to_insert_in_edge.range.start] as char, new_node_index + 1), // index is 1 higher than the already pushed top
                (new_node_char, new_node_index)
            ]
        );

        current_node.range.end = current_node.range.start + self.index;
        self.tree.arena.push(node_to_insert_in_edge); // execute the push!
    }

    pub fn add_leaf(&mut self, index_in_entry: usize, end_index: usize, input_string: &[u8]) {
        let new_node = Node::build_in_arena(Range::new(index_in_entry, end_index), self.current_node_index_in_arena, std::array::from_fn(|_| None), None, None);
        let new_node_index = self.tree.arena.len();
        self.tree.arena.push(new_node);
        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        current_node.add_child(input_string[index_in_entry] as char, new_node_index);
    }

}