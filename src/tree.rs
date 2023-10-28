use crate::tree_builder::TreeBuilder;
// TODO: kijk hoeveel nodes er gemaakt worden bij uit hetvoeren van de build (leaves en interne nodes, zoals in de cpp)
pub const MAX_CHILDREN: usize = 28;

/// Custom trait implemented by types that have a value that represents NULL
pub trait Nullable<T> {
    const NULL: T;

    fn is_null(&self) -> bool;
}

/// Type that represents the index of a node in the arena part of the tree
pub type NodeIndex = usize;

impl Nullable<NodeIndex> for NodeIndex {
    /// Use usize::MAX as NULL value since this will in practice never be reached.
    /// It is not possible to create 2^64-1 nodes (on a 64-bit machine). 
    /// This would simply never fit in memory
    const NULL: NodeIndex = usize::MAX;

    fn is_null(&self) -> bool {
        *self == NodeIndex::NULL
    }
}


#[derive(Debug, PartialEq)]
pub struct Tree {
    pub arena: Vec<Node>,
}

impl Tree {
    pub fn new(data: &str, builder: impl TreeBuilder) -> Self {
        builder.build(
            data,
            Tree {
                arena: vec![Node::create_root()],
            },
        )
    }
}

#[derive(Debug, PartialEq)]
pub struct Node {
    pub range: Range,
    pub children: [NodeIndex; MAX_CHILDREN],
    pub parent: NodeIndex,
    pub link: NodeIndex,
    pub suffix_index: NodeIndex,
}

impl Node {
    pub fn create_root() -> Self {
        Node {
            range: Range::new(0, 0),
            children: [NodeIndex::NULL; MAX_CHILDREN],
            parent: NodeIndex::NULL,
            link: NodeIndex::NULL,
            suffix_index: NodeIndex::NULL,
        }
    }

    /// Returns a tuple that contains the index of the new node in the arena and a reference to that node
    pub fn new(range: Range, parent: NodeIndex, children: [NodeIndex; MAX_CHILDREN], link: NodeIndex, suffix_index: usize) -> Node {
        Node {
            range,
            children,
            parent,
            link,
            suffix_index,
        }
    }

    pub fn new_with_child_tuples(range: Range, parent: NodeIndex, children_tuples: Vec<(u8, NodeIndex)>, link: NodeIndex, suffix_index: usize) -> Node {
        let mut node = Node {
            range,
            children: [NodeIndex::NULL; MAX_CHILDREN],
            parent,
            link,
            suffix_index,
        };
        children_tuples.iter().for_each(|(char, child)| node.add_child(*char, *child));
        node
    }

    fn char_to_child_index(character: u8) -> usize {
        if character == b'#' {
            26
        } else if character == b'$' {
            27
        } else {
            character as usize - 65 // 65 is 'A' in ascii
        }
    }

    pub fn add_child(&mut self, character: u8, child: NodeIndex) {
        self.children[Self::char_to_child_index(character)] = child;
    }

    // TODO: mogelijks rekening houden met feit dat er tekens ingeput kunnen worden die niet in de array zitten?
    //  (bv ?*^), dit gaat er nu voor zorgen dat alles crasht als we zo'n teken mee geven
    pub fn get_child(&self, character: u8) -> NodeIndex {
        self.children[Self::char_to_child_index(character)]
    }

    pub fn set_new_children(&mut self, new_children: Vec<(u8, NodeIndex)>) {
        self.children = [NodeIndex::NULL; MAX_CHILDREN];
        new_children.iter().for_each(|(character, child)| self.add_child(*character, *child));
    }
}

#[derive(Debug, PartialEq)]
pub struct Range {
    pub start: usize,
    pub end: usize,
}

impl Range {
    pub fn new(start: usize, end: usize) -> Self {
        Range { start, end }
    }
    pub fn length(&self) -> usize {
        self.end - self.start
    }
}

#[cfg(test)]
mod tests {
    use crate::tree::{MAX_CHILDREN, Node, NodeIndex, Nullable, Range, Tree};
    use crate::tree_builder::{TreeBuilder, UkkonenBuilder};

    #[test]
    fn total_test() {
        let input = "ACACACGT$";

        let tree = Tree::new(input, UkkonenBuilder::new());
        let mut control_tree = Tree { arena: vec![]};
        for _ in 0..14 {
            control_tree.arena.push(Node::new(
                Range::new(0, 0),
                NodeIndex::NULL,
                [NodeIndex::NULL; MAX_CHILDREN],
                NodeIndex::NULL,
                NodeIndex::NULL
            ));
        }

        // too see the required structure: place the input in: https://brenden.github.io/ukkonen-animation/

        // set the parents right
        control_tree.arena[1].parent = 3;
        control_tree.arena[2].parent = 5;
        control_tree.arena[3].parent = 7;
        control_tree.arena[4].parent = 3;
        control_tree.arena[5].parent = 9;
        control_tree.arena[6].parent = 5;
        control_tree.arena[7].parent = 0;
        control_tree.arena[8].parent = 7;
        control_tree.arena[9].parent = 0;
        control_tree.arena[10].parent = 9;
        control_tree.arena[11].parent = 0;
        control_tree.arena[12].parent = 0;
        control_tree.arena[13].parent = 0;

        // set children
        let children_for_each_node = vec![
            vec![('$', 13), ('A', 7), ('C', 9), ('G', 11), ('T', 12)],
            vec![],
            vec![],
            vec![('A', 1), ('G', 4)],
            vec![],
            vec![('A', 2), ('G', 6)],
            vec![],
            vec![('A', 3), ('G', 8)],
            vec![],
            vec![('A', 5), ('G', 10)],
            vec![],
            vec![],
            vec![],
            vec![],
        ];
        for (i, children) in children_for_each_node.iter().enumerate() {
            children.iter().for_each(|(character, child)| {
                control_tree.arena[i].add_child(*character as u8, *child);
            });
        }

        // set ranges
        control_tree.arena[1].range = Range::new(4, 9);
        control_tree.arena[2].range = Range::new(4, 9);
        control_tree.arena[3].range = Range::new(2, 4);
        control_tree.arena[4].range = Range::new(6, 9);
        control_tree.arena[5].range = Range::new(2, 4);
        control_tree.arena[6].range = Range::new(6, 9);
        control_tree.arena[7].range = Range::new(0, 2);
        control_tree.arena[8].range = Range::new(6, 9);
        control_tree.arena[9].range = Range::new(1, 2);
        control_tree.arena[10].range = Range::new(6, 9);
        control_tree.arena[11].range = Range::new(6, 9);
        control_tree.arena[12].range = Range::new(7, 9);
        control_tree.arena[13].range = Range::new(8, 9);

        // set suffix links
        control_tree.arena[3].link = 5;
        control_tree.arena[5].link = 7;
        control_tree.arena[7].link = 9;
        control_tree.arena[9].link = 0;

        // set suffix indices
        control_tree.arena[1].suffix_index = 0;
        control_tree.arena[4].suffix_index = 2;
        control_tree.arena[8].suffix_index = 4;
        control_tree.arena[2].suffix_index = 1;
        control_tree.arena[6].suffix_index = 3;
        control_tree.arena[10].suffix_index = 5;
        control_tree.arena[11].suffix_index = 6;
        control_tree.arena[12].suffix_index = 7;
        control_tree.arena[13].suffix_index = 8;

        assert_eq!(tree, control_tree);
    }
}