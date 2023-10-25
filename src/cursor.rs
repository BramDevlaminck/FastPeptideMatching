use crate::Node;
use crate::tree::Tree;

#[derive(Debug, PartialEq)]
pub enum CursorIterator {
    Ok,
    AtEnd,
    InWord
}

pub struct Cursor<'a> {
    pub current_node_index_in_arena: usize,
    pub current_node: &'a mut Node,
    pub index: usize,
    pub tree: &'a mut Tree
}

impl<'a> Cursor<'a> {

    pub fn new(tree: &'a mut Tree) -> Cursor<'a> {
        Cursor {
            current_node_index_in_arena: 0,
            current_node: &mut tree.arena[0],
            index: 0,
            tree
        }
    }

    pub fn next(&mut self, next_character: char) -> CursorIterator {
        if self.index < self.current_node.range.length() {
            if self.tree.string.as_bytes()[self.current_node.range.start] as char == next_character {
                self.index += 1;
                return CursorIterator::Ok;
            }
            return CursorIterator::InWord
        }

        if let Some(child) = self.current_node.get_child(next_character) {
            self.current_node_index_in_arena = child;
            self.current_node = &mut self.tree.arena[self.current_node_index_in_arena];
            self.index = 1;
            return CursorIterator::Ok
        }

        return CursorIterator::AtEnd
    }

    pub fn reset(&mut self) {
        self.index = 0;
        self.current_node = &mut self.tree.arena[0];
    }

}