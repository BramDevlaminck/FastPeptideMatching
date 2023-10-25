use std::cmp::min;
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

    /// the split function for the naive builder
    pub fn split_and_add_naive(&mut self, index_in_entry: usize, end_index: usize, input_string: &[u8]) {
        let new_node = Node::new(Range::new(index_in_entry, end_index), self.current_node_index_in_arena, std::array::from_fn(|_| None), None, None);
        let new_node_char = input_string[new_node.range.start] as char;
        let new_node_index = self.tree.arena.len();
        self.tree.arena.push(new_node);

        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        let node_to_insert_in_edge = Node::new(Range::new(current_node.range.start + self.index, current_node.range.end), self.current_node_index_in_arena, current_node.children, None, None);
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
        let new_node = Node::new(Range::new(index_in_entry, end_index), self.current_node_index_in_arena, std::array::from_fn(|_| None), None, None);
        let new_node_index = self.tree.arena.len();
        self.tree.arena.push(new_node);
        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        current_node.add_child(input_string[index_in_entry] as char, new_node_index);
    }

    pub fn at_node(&self) -> bool {
        self.index == self.tree.arena[self.current_node_index_in_arena].range.length()
    }

    pub fn add_link(&mut self, receiver: usize, link_to: usize) {
        self.tree.arena[receiver].link = Some(link_to);
    }

    /// Split edge implementation for Ukkonen
    pub fn split_edge(&mut self, input_string: &[u8]) -> usize {
        // first get the index where the next node will be inserted, do this before we have a mutable borrow
        let new_internal_node_index_in_arena = self.tree.arena.len();
        // create the new node
        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        let new_internal_node_end = current_node.range.start + self.index;
        let new_internal_node = Node::new_with_child_tuples(
            Range::new(current_node.range.start, new_internal_node_end),
            current_node.parent.unwrap(),
            vec![
                (input_string[new_internal_node_end] as char, self.current_node_index_in_arena)
            ],
            None,
            None
        );
        let parent_index_in_arena = current_node.parent.unwrap(); // temp store the index since we will need it later
        // update current node
        current_node.range.start += self.index;
        current_node.parent = Some(new_internal_node_index_in_arena);
        // update the parent now we have updated everything needed to the current node
        let parent = &mut self.tree.arena[parent_index_in_arena];
        parent.add_child(input_string[new_internal_node.range.start] as char, new_internal_node_index_in_arena);
        // actually push the new internal node and update the cursor
        self.tree.arena.push(new_internal_node);
        self.current_node_index_in_arena = new_internal_node_index_in_arena;

        new_internal_node_index_in_arena
    }

    /// Add a leaf with suffix index used in the Ukkonen implementation
    pub fn add_leaf_from_position(&mut self, j: usize, suffix_index: usize, input_string: &[u8]) {
        let new_leaf = Node::new(
            Range::new(j, input_string.len()),
            self.current_node_index_in_arena,
            std::array::from_fn(|_| None),
            None,
            Some(suffix_index)
        );
        let new_leaf_position_in_arena = self.tree.arena.len();
        let current_node = &mut self.tree.arena[self.current_node_index_in_arena];
        current_node.add_child(input_string[j] as char, new_leaf_position_in_arena);
        self.tree.arena.push(new_leaf);
    }

    pub fn follow_link(&mut self, input_string: &[u8]) {
        if self.current_node_index_in_arena == 0 {
            return;
        }

        let mut current_node = &self.tree.arena[self.current_node_index_in_arena];
        let mut begin = current_node.range.start;

        let mut distance_left_to_walk;
        if current_node.parent.unwrap() == 0 {
            self.current_node_index_in_arena = 0;
            begin += 1;
            distance_left_to_walk = self.index - 1;
        } else {
            // follow link
            distance_left_to_walk = self.index; // distance before following link
            self.current_node_index_in_arena = self.tree.arena[current_node.parent.unwrap()].link.unwrap();
        }
        current_node = &self.tree.arena[self.current_node_index_in_arena];
        self.index = current_node.range.length();

        while distance_left_to_walk > 0 {
            // move to child
            self.current_node_index_in_arena = current_node.get_child(input_string[begin] as char).unwrap();
            current_node = &self.tree.arena[self.current_node_index_in_arena];

            // walk as far as possible on current edge
            let current_advance = min(current_node.range.length(), distance_left_to_walk);
            distance_left_to_walk -= current_advance;
            begin += current_advance;
            self.index = current_advance;
        }
    }

}