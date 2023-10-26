use crate::cursor::CursorIterator;
use crate::tree::{NULL, Tree};

pub struct ReadOnlyCursor<'a> {
    pub current_node_index_in_arena: usize,
    pub index: usize,
    pub tree: &'a Tree
}

impl<'a> ReadOnlyCursor<'a> {
    pub fn new(tree: &'a Tree) -> ReadOnlyCursor<'a> {
        Self {
            current_node_index_in_arena: 0,
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

        let child = current_node.get_child(next_character);
        if child != NULL {
            self.current_node_index_in_arena = child;
            self.index = 1;
            return CursorIterator::Ok
        }

        CursorIterator::AtEnd
    }

    pub fn reset(&mut self) {
        self.index = 0;
        self.current_node_index_in_arena = 0;
    }


}