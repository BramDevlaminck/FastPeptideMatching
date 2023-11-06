use umgap::{agg, rmq::lca::LCACalculator, taxon};
use umgap::agg::Aggregator;
use umgap::taxon::{TaxonId, TaxonList, TaxonTree};

use crate::Protein;
use crate::tree::{NodeIndex, Nullable, Tree};


static mut PRINT_INDEX: i32 = 1;


pub struct TaxonIdCalculator {
    snapping: Vec<Option<TaxonId>>,
    aggregator: Box<dyn Aggregator>,
}

impl TaxonIdCalculator {
    pub fn new(ncbi_taxonomy_fasta_file: &str) -> Box<Self> {
        let taxons = taxon::read_taxa_file(ncbi_taxonomy_fasta_file).unwrap();
        let taxon_tree = TaxonTree::new(&taxons);
        let by_id = TaxonList::new(taxons);
        let snapping = taxon_tree.snapping(&by_id, true);

        let aggregator = LCACalculator::new(taxon_tree);

        Box::new(Self {
            snapping,
            aggregator: Box::new(aggregator),
        })
    }

    /// Calculates the taxon ids by only using the leaves in the tree
    pub fn calculate_taxon_ids_leaf(&self, tree: &mut Tree, proteins: &Vec<Protein>) {
        self.calculate_taxon_ids_leaf_recursive(tree, proteins, 0);
    }

    /// Calculate the taxon id for current_node by only using the leaves in the tree
    fn calculate_taxon_ids_leaf_recursive(&self, tree: &mut Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) {
        let current_node = &mut tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            current_node.taxon_id = proteins[current_node.suffix_index].id;
            return
        }

        let taxon_id = self.get_aggregate(Self::get_taxon_id_leaves_under_node(tree, proteins, current_node_index));

        let current_node = &mut tree.arena[current_node_index];
        current_node.taxon_id = taxon_id;

        for child in current_node.children {
            if !child.is_null() {
                self.calculate_taxon_ids_leaf_recursive(tree, proteins, child);
            }
        }
    }
    /// Fill in all the taxon ids for the complete tree
    pub fn calculate_taxon_ids(&self, tree: &mut Tree, proteins: &Vec<Protein>) {
        self.recursive_calculate_taxon_ids(tree, proteins, 0);
    }

    /// The recursive function called by `calculate_taxon_ids`
    /// This function traverses the tree post order and fills in the taxon ids
    fn recursive_calculate_taxon_ids(&self, tree: &mut Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) -> TaxonId {
        let current_node = &mut tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            current_node.taxon_id = proteins[current_node.suffix_index].id;
            return current_node.taxon_id;
        }

        let mut taxon_ids = vec![];
        for child in current_node.children {
            if !child.is_null() {
                taxon_ids.push(self.recursive_calculate_taxon_ids(tree, proteins, child));
            }
        }

        let taxon_id = self.get_aggregate(taxon_ids);

        let current_node = &mut tree.arena[current_node_index];
        current_node.taxon_id = taxon_id;
        taxon_id
    }

    /// returns the taxon ids of all the leaves that are under current_node
    fn get_taxon_id_leaves_under_node(tree: &Tree, proteins: &Vec<Protein>, current_node_index: NodeIndex) -> Vec<TaxonId>{
        let mut ids: Vec<TaxonId> = vec![];
        Self::get_taxon_id_leaves_recursive(tree, proteins, &mut ids, current_node_index);
        ids
    }

    /// recursive function that adds al the taxon ids of the leaves under current_node to taxon_ids
    fn get_taxon_id_leaves_recursive(tree: &Tree, proteins: &Vec<Protein>, taxon_ids: &mut Vec<TaxonId>, current_node_index: NodeIndex) {
        let current_node = &tree.arena[current_node_index];
        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            taxon_ids.push(proteins[current_node.suffix_index].id);
            return
        }

        for child in current_node.children {
            if !child.is_null() {
                Self::get_taxon_id_leaves_recursive(tree, proteins, taxon_ids, child);
            }
        }
    }

    /// Helper function to print out the taxon ids of the tree in pre-order traversal
    pub fn output_all_taxon_ids(tree: &Tree, proteins: &Vec<Protein>) {
        Self::output_all_taxon_ids_recursive(tree, 0, proteins)
    }

    fn output_all_taxon_ids_recursive(tree: &Tree, current_node_index: TaxonId, proteins: &Vec<Protein>) {
        let current_node = &tree.arena[current_node_index];

        // we are in a leaf
        if !current_node.suffix_index.is_null() {
            println!("{}", current_node.taxon_id);
            println!("children: /");
            unsafe {
                PRINT_INDEX += 1;
            }
            return
        }

        for child in current_node.children {
            if !child.is_null() {
                Self::output_all_taxon_ids_recursive(tree, child, proteins);
            }
        }

        println!("{}", current_node.taxon_id);
        unsafe {
            if PRINT_INDEX == 21085 {
                let mut taxon_ids = vec![];
                Self::get_taxon_id_leaves_recursive(tree, proteins, &mut taxon_ids, current_node_index);
                println!("leaves: {:?}", taxon_ids);
            }
        }
        print!("children: ");
        for child in current_node.children {
            if !child.is_null() {
                print!("{};", tree.arena[child].taxon_id);
            }
        }
        println!();
        unsafe {
            PRINT_INDEX += 1;
        }
    }

    /// Snaps the given taxon ID using the ncbi taxonomy that was provided to the TaxonIdCalculator.
    pub fn snap_taxon_id(&self, id: TaxonId) -> TaxonId {
        self.snapping[id].unwrap_or_else(|| panic!("Could not snap taxon with id {id}"))
    }

    /// Calculates the aggregate of a vector containing TaxonIds
    fn get_aggregate(&self, ids: Vec<TaxonId>) -> TaxonId {
        let count = agg::count(ids.into_iter().filter(|&id| id != 0).map(|it| (it, 1.0))); // TODO: waarom de filter id != 0 ?
        self.aggregator.aggregate(&count).unwrap_or_else(|_| panic!("Could not aggregate following taxon ids: {:?}", &count))
    }
}


#[cfg(test)]
mod test {
    use crate::taxon_id_calculator::TaxonIdCalculator;
    use crate::Protein;
    use crate::tree::{MAX_CHILDREN, Node, NodeIndex, Nullable, Range, Tree};

    #[test]
    fn test_calculate_taxon_ids() {
        // the tree structure we are building in this test
        // with the expected taxon id between parentheses under the node id
        //              0
        //             (1)
        //             /|\
        //            / | \
        //           /  |  \
        //          /   |   \
        //         1    2    3
        //        (6)  (20) (2)
        //        / \
        //       /   \
        //      4     5
        //     (10)  (9)

        let test_taxonomy_file = "small_taxonomy.tsv";
        let mut tree = Tree {
            arena: vec![
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
                Node {
                    range: Range::new(0, 0),
                    children: [NodeIndex::NULL; MAX_CHILDREN],
                    parent: NodeIndex::NULL,
                    link: NodeIndex::NULL,
                    suffix_index: NodeIndex::NULL,
                    taxon_id: 1,
                },
            ]
        };

        // set some child structure
        tree.arena[0].children[0] = 1;
        tree.arena[0].children[3] = 2;
        tree.arena[0].children[5] = 3;
        tree.arena[1].children[1] = 4;
        tree.arena[1].children[5] = 5;

        // set the parents to match what we set in the children
        tree.arena[1].parent = 0;
        tree.arena[2].parent = 0;
        tree.arena[3].parent = 0;
        tree.arena[4].parent = 1;
        tree.arena[5].parent = 1;

        // set suffix ids in the leaves
        tree.arena[4].suffix_index = 0;
        tree.arena[5].suffix_index = 1;
        tree.arena[2].suffix_index = 2;
        tree.arena[3].suffix_index = 3;


        let proteins = vec![
            Protein {
                sequence: "ABC".to_string(),
                id: 10,
            },
            Protein {
                sequence: "XYZ".to_string(),
                id: 9,
            },
            Protein {
                sequence: "XAZ".to_string(),
                id: 20,
            },
            Protein {
                sequence: "W".to_string(),
                id: 2,
            },
        ];

        TaxonIdCalculator::new(test_taxonomy_file).calculate_taxon_ids(&mut tree, &proteins);

        assert_eq!(tree.arena[0].taxon_id, 1);
        assert_eq!(tree.arena[1].taxon_id, 6);
        assert_eq!(tree.arena[2].taxon_id, 20);
        assert_eq!(tree.arena[3].taxon_id, 2);
        assert_eq!(tree.arena[4].taxon_id, 10);
        assert_eq!(tree.arena[5].taxon_id, 9);
    }
}