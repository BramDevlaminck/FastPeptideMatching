use crate::tree_builder::TreeBuilder;

pub const NULL: usize = usize::MAX;

pub const MAX_CHILDREN: usize = 28;

#[derive(Debug)]
pub struct Tree {
    pub arena: Vec<Node>,
}

impl Tree {
    pub fn new(data: &String, builder: impl TreeBuilder) -> Self {
        builder.build(
            data,
            Tree {
                arena: vec![Node::create_root()],
            }
        )
    }
}

#[derive(Debug)]
pub struct Node {
    pub range: Range,
    pub children: [usize; MAX_CHILDREN],
    pub parent: usize,
    pub link: usize,
    pub suffix_index: usize,
}

impl Node {

    pub fn create_root() -> Self {
        Node {
            range: Range::new(0, 0),
            children: [NULL; MAX_CHILDREN],
            parent: NULL,
            link: NULL,
            suffix_index: NULL,
        }
    }

    /// Returns a tuple that contains the index of the new node in the arena and a reference to that node
    pub fn new(range: Range, parent: usize, children: [usize; MAX_CHILDREN], link: usize, suffix_index: usize) -> Node {
        Node {
            range,
            children,
            parent,
            link,
            suffix_index,
        }
    }

    pub fn new_with_child_tuples(range: Range, parent: usize, children_tuples: Vec<(char, usize)>, link: usize, suffix_index: usize) -> Node {
        let mut node = Node {
            range,
            children: [NULL; MAX_CHILDREN],
            parent,
            link,
            suffix_index,
        };
        children_tuples.iter().for_each(|(char, child)| node.add_child(*char, *child));
        node
    }

    fn char_to_child_index(character: char) -> usize {
        if character == '#' {
            26
        } else if character == '$' {
            27
        } else {
            character as usize - 'A' as usize
        }
    }

    pub fn add_child(&mut self, character: char, child: usize) {
        self.children[Self::char_to_child_index(character)] = child;
    }

    // TODO: mogelijks rekening houden met feit dat er tekens ingeput kunnen worden die niet in de array zitten?
    //  (bv ?*^), dit gaat er nu voor zorgen dat alles crasht als we zo'n teken mee geven
    pub fn get_child(&self, character: char) -> usize {
        self.children[Self::char_to_child_index(character)]
    }

    pub fn set_new_children(&mut self, new_children: Vec<(char, usize)>) {
        self.children = [NULL; MAX_CHILDREN];
        new_children.iter().for_each(|(character, child)| self.add_child(*character, *child));
    }
}

#[derive(Debug)]
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