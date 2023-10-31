use umgap::{agg, taxon, tree::lca::LCACalculator};
use umgap::agg::Aggregator;
use umgap::taxon::TaxonId;

use crate::Protein;
use crate::tree::{NodeIndex, Nullable, Tree};

pub fn calculate_lcas(tree: &mut Tree, ncbi_taxonomy_fasta_file: &String, proteins: &Vec<Protein>) {
    let taxons = taxon::read_taxa_file(ncbi_taxonomy_fasta_file).unwrap(); // TODO: handle failure
    let taxon_tree = taxon::TaxonTree::new(&taxons);
    let by_id = taxon::TaxonList::new(taxons);
    let snapping = taxon_tree.snapping(&by_id, false);

    let calculator = LCACalculator::new(taxon_tree.root, &by_id);

    recursive_calculate_lcas(tree, proteins, &calculator, &snapping, 0);
}

fn recursive_calculate_lcas(tree: &mut Tree, proteins: &Vec<Protein>, calculator: &LCACalculator, snapping: &Vec<Option<TaxonId>>, current_node_index: NodeIndex) -> TaxonId {
    let current_node = &mut tree.arena[current_node_index];
    // we are in a leave
    if !current_node.suffix_index.is_null() {
        current_node.taxon_id = proteins[current_node.suffix_index].id;
        return current_node.taxon_id;
    }

    let mut taxon_ids = vec![];
    for child in current_node.children {
        if !child.is_null() {
            taxon_ids.push(recursive_calculate_lcas(tree, proteins, calculator, snapping, child));
        }
    }

    let taxon_id = get_lca(calculator, snapping, taxon_ids);

    let current_node = &mut tree.arena[current_node_index];
    current_node.taxon_id = taxon_id;
    taxon_id
}


fn get_lca(calculator: &LCACalculator, snapping: &[Option<TaxonId>], ids: Vec<TaxonId>) -> TaxonId {
    let count = agg::count(ids.into_iter().map(|it| (it, 1.0)).filter(|(id, _)| *id != 0));
    // let counts = agg::filter(counts, args.lower_bound); TODO: used in umgap, but probably not needed here?
    let aggregate = calculator.aggregate(&count).unwrap(); // TODO: handle errors
    snapping[aggregate].unwrap()
}