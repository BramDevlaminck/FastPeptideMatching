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

    pub fn find_in_tree(&mut self, tree_data: &[u8], search_data: &[u8]) -> (bool, NodeIndex) {
        if search_data.is_empty() {
            return (true, 0);
        }

        let search_data_len = search_data.len();
        let mut current_node = &self.tree.arena[self.current_node_index_in_arena];
        let mut child_to_move_to = current_node.get_child(search_data[0]); // get us started on the right path
        let mut current_index_in_search_data = 0;
        while !child_to_move_to.is_null() && current_index_in_search_data <= search_data_len {
            current_node = &self.tree.arena[child_to_move_to];
            let min_dist = min(current_node.range.length(), search_data_len - current_index_in_search_data);
            // if there is a difference somewhere in the node => no match
            let next_index_in_search_data = current_index_in_search_data + min_dist;
            if search_data[current_index_in_search_data..next_index_in_search_data] != tree_data[current_node.range.start..current_node.range.start+min_dist] {
                return (false, child_to_move_to);
            }

            // at and of search string => we have a match
            current_index_in_search_data = next_index_in_search_data;
            if current_index_in_search_data == search_data_len {
                return (true, child_to_move_to);
            }

            // move to next child
            child_to_move_to = current_node.get_child(search_data[current_index_in_search_data]);
        }

        (false, child_to_move_to)
    }

    pub fn reset(&mut self) {
        self.index = 0;
        self.current_node_index_in_arena = 0;
    }
}