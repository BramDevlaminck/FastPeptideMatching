use std::cmp::min;
use crate::cursor::CursorIterator;
use crate::tree::{NodeIndex, Nullable, Tree};

/// A Cursor that cannot mutate the tree (which means it can only be used during the search phase)
pub struct ReadOnlyCursor<'a> {
    pub current_node_index_in_arena: usize,
    pub index: usize,
    pub tree: &'a Tree,
}

impl<'a> ReadOnlyCursor<'a> {
    pub fn new(tree: &'a Tree) -> ReadOnlyCursor<'a> {
        Self {
            current_node_index_in_arena: 0,
            index: 0,
            tree,
        }
    }

    /// Try to progress by consuming `next_character`
    /// Returns Some(()) if we were able to move to the next location.
    /// None otherwise
    pub fn next(&mut self, next_character: u8, bytes_input: &[u8]) -> Option<()> {
        let current_node = &self.tree.arena[self.current_node_index_in_arena];
        if self.index < current_node.range.length() {
            if bytes_input[current_node.range.start + self.index] == next_character {
                self.index += 1;
                return Some(());
            }
            return None;
        }

        let child = current_node.get_child(next_character);
        if !child.is_null() {
            self.current_node_index_in_arena = child;
            self.index = 1;
            return Some(());
        }

        None
    }

    pub fn reset(&mut self) {
        self.index = 0;
        self.current_node_index_in_arena = 0;
    }
}