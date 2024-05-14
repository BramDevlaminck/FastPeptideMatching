use clap::ValueEnum;
use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
use crate::Nullable;

/// Enum used to define the commandline arguments and choose which index style is used
#[derive(ValueEnum, Clone, Debug, PartialEq)]
pub enum SuffixToProteinMappingStyle {
    Dense,
    Sparse
}

/// Trait implemented by the SuffixToProtein mappings
pub trait SuffixToProteinIndex: Send + Sync {
    
    /// Returns the index of the protein in the protein list for the given suffix
    ///
    /// # Arguments
    /// * `suffix` - The suffix of which we want to know of which protein it is a part
    ///
    /// # Returns
    ///
    /// Returns the index of the protein in the proteins list of which the suffix is a part
    fn suffix_to_protein(&self, suffix: i64) -> u32;
}

/// Mapping that uses O(n) memory with n the size of the input text, but retrieval of the protein is in O(1)
#[derive(Debug, PartialEq)]
pub struct DenseSuffixToProtein {
    // UniProtKB does not have more that u32::MAX proteins, so a larger type is not needed
    mapping: Vec<u32>,
}

/// Mapping that uses O(m) memory with m the number of proteins, but retrieval of the protein is O(log m)
#[derive(Debug, PartialEq)]
pub struct SparseSuffixToProtein {
    mapping: Vec<i64>,
}

impl SuffixToProteinIndex for DenseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        self.mapping[suffix as usize]
    }
}

impl SuffixToProteinIndex for SparseSuffixToProtein {
    fn suffix_to_protein(&self, suffix: i64) -> u32 {
        let protein_index = self.mapping.binary_search(&suffix).unwrap_or_else(|index| index - 1);
        // if the next value in the mapping is 1 larger than the current suffix, that means that the current suffix starts with a SEPARATION_CHARACTER or TERMINATION_CHARACTER
        // this means it does not belong to a protein
        if self.mapping[protein_index + 1] == suffix + 1 {
            return u32::NULL
        }
        protein_index as u32
    }
}

impl DenseSuffixToProtein {

    /// Creates a new DenseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new DenseSuffixToProtein build over the provided text
    pub fn new(text: &[u8]) -> Self {
        let mut current_protein_index: u32 = 0;
        let mut suffix_index_to_protein: Vec<u32> = vec![];
        for &char in text.iter() {
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                current_protein_index += 1;
                suffix_index_to_protein.push(u32::NULL);
            } else {
                assert_ne!(current_protein_index, u32::NULL);
                suffix_index_to_protein.push(current_protein_index);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        DenseSuffixToProtein { mapping: suffix_index_to_protein }
    }
}

impl SparseSuffixToProtein {

    /// Creates a new SparseSuffixToProtein mapping
    ///
    /// # Arguments
    /// * `text` - The text over which we want to create the mapping
    ///
    /// # Returns
    ///
    /// Returns a new SparseSuffixToProtein build over the provided text
    pub fn new(text: &[u8]) -> Self {
        let mut suffix_index_to_protein: Vec<i64> = vec![0];
        for (index, &char) in text.iter().enumerate() {
            if char == SEPARATION_CHARACTER || char == TERMINATION_CHARACTER {
                suffix_index_to_protein.push(index as i64 + 1);
            }
        }
        suffix_index_to_protein.shrink_to_fit();
        SparseSuffixToProtein { mapping: suffix_index_to_protein }
    }

}


#[cfg(test)]
mod tests {
    use sa_mappings::proteins::{SEPARATION_CHARACTER, TERMINATION_CHARACTER};
    use crate::Nullable;
    use crate::suffix_to_protein_index::{DenseSuffixToProtein, SparseSuffixToProtein, SuffixToProteinIndex};

    fn build_text() -> Vec<u8> {
        let mut text = ["ACG", "CG", "AAA"].join(&format!("{}", SEPARATION_CHARACTER as char));
        text.push(TERMINATION_CHARACTER as char);
        text.into_bytes()
    }

    #[test]
    fn test_dense_build() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        let expected = DenseSuffixToProtein {mapping: vec![0, 0, 0, u32::NULL, 1, 1, u32::NULL, 2, 2, 2, u32::NULL]};
        assert_eq!(index, expected);
    }

    #[test]
    fn test_sparse_build() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        let expected = SparseSuffixToProtein {mapping: vec![0, 4, 7, 11]};
        assert_eq!(index, expected);
    }

    #[test]
    fn test_search_dense() {
        let u8_text = &build_text();
        let index = DenseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }

    #[test]
    fn test_search_sparse() {
        let u8_text = &build_text();
        let index = SparseSuffixToProtein::new(u8_text);
        assert_eq!(index.suffix_to_protein(5), 1);
        assert_eq!(index.suffix_to_protein(7), 2);
        // suffix that starts with SEPARATION_CHARACTER
        assert_eq!(index.suffix_to_protein(3), u32::NULL);
        // suffix that starts with TERMINATION_CHARACTER
        assert_eq!(index.suffix_to_protein(10), u32::NULL);
    }
}