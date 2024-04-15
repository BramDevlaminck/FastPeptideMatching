//! This module provides a `TaxonAggregator` struct that is used to aggregate taxonomic information.
//! It uses a taxonomy file to create a taxonomic tree and performs aggregation using different
//! methods.

use std::error::Error;

use umgap::{
    agg::{
        count,
        MultiThreadSafeAggregator
    },
    rmq::{
        lca::LCACalculator,
        mix::MixCalculator
    },
    taxon::{
        read_taxa_file,
        TaxonId,
        TaxonList,
        TaxonTree
    }
};

/// A struct that represents a taxon aggregator.
pub struct TaxonAggregator {
    /// A vector that contains the snapped taxon IDs.
    snapping: Vec<Option<TaxonId>>,

    /// The aggregator used to aggregate taxon IDs.
    aggregator: Box<dyn MultiThreadSafeAggregator>,

    /// The taxon list.
    taxon_list: TaxonList
}

/// An enum that specifies the aggregation method to use.
pub enum AggregationMethod {
    /// The Lowest Common Ancestor (LCA) aggregation method.
    Lca,

    /// The LCA* aggregation method.
    LcaStar
}

impl TaxonAggregator {
    /// Creates a new `TaxonAggregator` from a taxonomy file and an aggregation method.
    ///
    /// # Arguments
    ///
    /// * `file` - A string slice that represents the path to the taxonomy file.
    /// * `method` - An `AggregationMethod` enum that specifies the aggregation method to use.
    ///
    /// # Returns
    ///
    /// Returns a `Result` containing the `TaxonAggregator`
    ///
    /// # Errors
    ///
    /// Returns a `Box<dyn Error>` if an error occurred while reading the taxonomy file.
    pub fn try_from_taxonomy_file(
        file: &str,
        method: AggregationMethod
    ) -> Result<Self, Box<dyn Error>> {
        let taxons = read_taxa_file(file)?;
        let taxon_tree = TaxonTree::new(&taxons);
        let taxon_list = TaxonList::new(taxons);
        let snapping = taxon_tree.snapping(&taxon_list, true);

        let aggregator: Box<dyn MultiThreadSafeAggregator> = match method {
            AggregationMethod::Lca => Box::new(MixCalculator::new(taxon_tree, 1.0)),
            AggregationMethod::LcaStar => Box::new(LCACalculator::new(taxon_tree))
        };

        Ok(Self {
            snapping,
            aggregator,
            taxon_list
        })
    }

    /// Checks if a taxon exists in the taxon list.
    ///
    /// # Arguments
    ///
    /// * `taxon` - The taxon ID to check.
    ///
    /// # Returns
    ///
    /// Returns a boolean value indicating whether the taxon exists in the taxon list.
    pub fn taxon_exists(&self, taxon: TaxonId) -> bool {
        self.taxon_list.get(taxon).is_some()
    }

    /// Snaps a taxon to its closest ancestor in the taxonomic tree.
    ///
    /// # Arguments
    ///
    /// * `taxon` - The taxon ID to snap.
    ///
    /// # Returns
    ///
    /// Returns the snapped taxon ID, or panics if the taxon cannot be snapped.
    pub fn snap_taxon(&self, taxon: TaxonId) -> TaxonId {
        self.snapping[taxon].unwrap_or_else(|| panic!("Could not snap taxon with id {taxon}"))
    }

    /// Aggregates a list of taxon IDs using the specified aggregation method.
    ///
    /// # Arguments
    ///
    /// * `taxa` - A vector of taxon IDs to aggregate.
    ///
    /// # Returns
    ///
    /// Returns the aggregated taxon ID, or panics if aggregation fails.
    pub fn aggregate(&self, taxa: Vec<TaxonId>) -> TaxonId {
        let count = count(taxa.into_iter().map(|t| (t, 1.0)));
        self.aggregator
            .aggregate(&count)
            .unwrap_or_else(|_| panic!("Could not aggregate following taxon ids: {:?}", &count))
    }
}

#[cfg(test)]
mod tests {
    use std::{
        fs::File,
        io::Write,
        path::PathBuf
    };

    use tempdir::TempDir;

    use super::*;

    fn create_taxonomy_file(tmp_dir: &TempDir) -> PathBuf {
        let taxonomy_file = tmp_dir.path().join("taxonomy.tsv");
        let mut file = File::create(&taxonomy_file).unwrap();

        writeln!(file, "1\troot\tno rank\t1\t\x01").unwrap();
        writeln!(file, "2\tBacteria\tsuperkingdom\t1\t\x01").unwrap();
        writeln!(file, "6\tAzorhizobium\tgenus\t1\t\x01").unwrap();
        writeln!(file, "7\tAzorhizobium caulinodans\tspecies\t6\t\x01").unwrap();
        writeln!(file, "9\tBuchnera aphidicola\tspecies\t6\t\x01").unwrap();
        writeln!(file, "10\tCellvibrio\tgenus\t6\t\x01").unwrap();
        writeln!(file, "11\tCellulomonas gilvus\tspecies\t10\t\x01").unwrap();
        writeln!(file, "13\tDictyoglomus\tgenus\t11\t\x01").unwrap();
        writeln!(file, "14\tDictyoglomus thermophilum\tspecies\t10\t\x01").unwrap();
        writeln!(file, "16\tMethylophilus\tgenus\t14\t\x01").unwrap();
        writeln!(file, "17\tMethylophilus methylotrophus\tspecies\t16\t\x01").unwrap();
        writeln!(file, "18\tPelobacter\tgenus\t17\t\x01").unwrap();
        writeln!(file, "19\tSyntrophotalea carbinolica\tspecies\t17\t\x01").unwrap();
        writeln!(file, "20\tPhenylobacterium\tgenus\t19\t\x01").unwrap();

        taxonomy_file
    }

    #[test]
    fn test_try_from_taxonomy_file() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_try_from_taxonomy_file").unwrap();

        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let _ = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();
        let _ = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::LcaStar
        )
        .unwrap();
    }

    #[test]
    fn test_taxon_exists() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_taxon_exists").unwrap();

        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();

        for i in 0 ..= 20 {
            if [0, 3, 4, 5, 8, 12, 15].contains(&i) {
                assert!(!taxon_aggregator.taxon_exists(i));
            } else {
                assert!(taxon_aggregator.taxon_exists(i));
            }
        }
    }

    #[test]
    fn test_snap_taxon() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_snap_taxon").unwrap();

        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();

        for i in 0 ..= 20 {
            if ![0, 3, 4, 5, 8, 12, 15].contains(&i) {
                assert_eq!(taxon_aggregator.snap_taxon(i), i);
            }
        }
    }

    #[test]
    fn test_aggregate_lca() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_aggregate").unwrap();

        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::Lca
        )
        .unwrap();

        assert_eq!(taxon_aggregator.aggregate(vec![7, 9]), 6);
        assert_eq!(taxon_aggregator.aggregate(vec![11, 14]), 10);
        assert_eq!(taxon_aggregator.aggregate(vec![17, 19]), 17);
    }

    #[test]
    fn test_aggregate_lca_star() {
        // Create a temporary directory for this test
        let tmp_dir = TempDir::new("test_aggregate").unwrap();

        let taxonomy_file = create_taxonomy_file(&tmp_dir);

        let taxon_aggregator = TaxonAggregator::try_from_taxonomy_file(
            taxonomy_file.to_str().unwrap(),
            AggregationMethod::LcaStar
        )
        .unwrap();

        assert_eq!(taxon_aggregator.aggregate(vec![7, 9]), 6);
        assert_eq!(taxon_aggregator.aggregate(vec![11, 14]), 10);
        assert_eq!(taxon_aggregator.aggregate(vec![17, 19]), 19);
    }
}
