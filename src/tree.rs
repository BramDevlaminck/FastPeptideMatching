use crate::tree_builder::TreeBuilder;

const MAX_CHILDREN: usize = 28;

#[derive(Debug)]
pub struct Tree {
    pub arena: Vec<Node>,
    // pub string: &'a String,
}

impl<'a> Tree {
    pub fn new(data: String, builder: impl TreeBuilder) -> Self {
        builder.build(
            &data,
            Tree {
                arena: vec![Node::create_root()],
            }
        )
    }
}

#[derive(Debug)]
pub struct Node {
    pub range: Range,
    pub children: [Option<usize>; MAX_CHILDREN],
    pub parent: Option<usize>,
    pub link: Option<usize>,
    pub suffix_index: Option<usize>,
}

impl Node {

    pub fn create_root() -> Self {
        Node {
            range: Range::new(0, 0),
            children: std::array::from_fn(|_| None),
            parent: None,
            link: None,
            suffix_index: None,
        }
    }

    /// Returns a tuple that contains the index of the new node in the arena and a reference to that node
    pub fn new<'a>(range: Range, parent: usize, children: [Option<usize>; MAX_CHILDREN], link: Option<usize>, suffix_index: Option<usize>) -> Node {
        Node {
            range,
            children,
            parent: Some(parent),
            link,
            suffix_index,
        }
    }
    
    pub fn new_with_child_tuples<'a>(range: Range, parent: usize, children_tuples: Vec<(char, usize)>, link: Option<usize>, suffix_index: Option<usize>) -> Node {
        let mut node = Node {
            range,
            children: std::array::from_fn(|_| None),
            parent: Some(parent),
            link,
            suffix_index,
        };
        children_tuples.iter().for_each(|(char, child)| node.add_child(*char, *child));
        node
    }

/*    pub fn new_boxed(range: Range, parent: i32, link: Option<i32>, suffix_index: Option<i32>) -> Box<Self> {
        Box::new(Self::new(range, parent, link, suffix_index))
    }*/

    fn char_to_child_index(character: char) -> usize {
        if character == '#' {
            26
        } else if character == '$' {
            27
        } else {
            character as usize - 'A' as usize
        }
    }

    // pub fn has_child(&self, character: char) -> bool {
    //     self.children.get(self.char_to_child_index(character)).unwrap().is_some()
    // }

    pub fn add_child(&mut self, character: char, child: usize) {
        self.children[Self::char_to_child_index(character)] = Some(child);
    }

    pub fn get_child(&self, character: char) -> Option<usize> {
        self.children[Self::char_to_child_index(character)]
    }

    pub fn set_new_children(&mut self, new_children: Vec<(char, usize)>) {
        self.children = std::array::from_fn(|_| None);
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