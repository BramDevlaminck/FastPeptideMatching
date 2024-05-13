use std::cmp::min;


use sa_mappings::functionality::{FunctionAggregator, FunctionalAggregation};
use sa_mappings::proteins::{Protein, Proteins};
use sa_mappings::taxonomy::TaxonAggregator;
use umgap::taxon::TaxonId;

use crate::sa_searcher::BoundSearch::{Maximum, Minimum};
use crate::suffix_to_protein_index::SuffixToProteinIndex;
use crate::Nullable;

/// Enum indicating if we are searching for the minimum, or maximum bound in the suffix array
#[derive(Clone, Copy, PartialEq)]
enum BoundSearch {
    Minimum,
    Maximum,
}

/// Enum representing the minimum and maximum bound of the found matches in the suffix array
#[derive(PartialEq, Debug)]
pub enum BoundSearchResult {
    NoMatches,
    SearchResult((usize, usize)),
}

/// Enum representing the matching suffixes after searching a peptide in the suffix array
/// Both the MaxMatches and SearchResult indicate found suffixes, but MaxMatches is used when the cutoff is reached.
#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    NoMatches,
    MaxMatches(Vec<i64>),
    SearchResult(Vec<i64>),
}

/// Custom implementation of partialEq for SearchAllSuffixesResult
/// We consider 2 SearchAllSuffixesResult equal if they exist of the same key, and the Vec contains the same values, but the order can be different
impl PartialEq for SearchAllSuffixesResult {
    fn eq(&self, other: &Self) -> bool {

        /// Returns true if `arr1` and `arr2` contains the same elements, the order of the elements is ignored
        ///
        /// # Arguments
        /// * `arr1` - The first array used in the comparison
        /// * `arr2` - The second array used in the comparison
        ///
        /// # Returns
        ///
        /// Returns true if arr1 and arr2 contains the same elements, the order of the elements is ignored
        fn array_eq_unordered(arr1: &[i64], arr2: &[i64]) -> bool {
            let mut arr1_copy = arr1.to_owned();
            let mut arr2_copy = arr2.to_owned();

            arr1_copy.sort();
            arr2_copy.sort();

            arr1_copy == arr2_copy
        }

        match (self, other) {
            (
                SearchAllSuffixesResult::MaxMatches(arr1),
                SearchAllSuffixesResult::MaxMatches(arr2),
            ) => array_eq_unordered(arr1, arr2),
            (
                SearchAllSuffixesResult::SearchResult(arr1),
                SearchAllSuffixesResult::SearchResult(arr2),
            ) => array_eq_unordered(arr1, arr2),
            (SearchAllSuffixesResult::NoMatches, SearchAllSuffixesResult::NoMatches) => true,
            _ => false,
        }
    }
}

/// Struct that contains all the elements needed to search a peptide in the suffix array
/// This struct also contains all the functions used for search
///
/// # Arguments
/// * `sa` - The sparse suffix array representing the protein database
/// * `sparseness_factor` - The sparseness factor used by the suffix array
/// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix is part of
/// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the taxonomic analysis provided by Unipept
/// * `function_aggregator` - Object used to retrieve the functional annotations and to calculate the functional analysis provided by Unipept
pub struct Searcher {
    sa: Vec<i64>,
    pub sparseness_factor: u8,
    suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
    proteins: Proteins,
    taxon_id_calculator: TaxonAggregator,
    function_aggregator: FunctionAggregator
}

impl Searcher {
    
    /// Creates a new Searcher object
    ///
    /// # Arguments
    /// * `sa` - The sparse suffix array representing the protein database
    /// * `sparseness_factor` - The sparseness factor used by the suffix array
    /// * `suffix_index_to_protein` - Mapping from a suffix to the proteins to know which a suffix is part of
    /// * `proteins` - List of all the proteins where the suffix array is build on
    /// * `taxon_id_calculator` - Object representing the used taxonomy and that calculates the taxonomic analysis provided by Unipept
    /// * `function_aggregator` - Object used to retrieve the functional annotations and to calculate the functional analysis provided by Unipept
    ///
    /// # Returns
    ///
    /// Returns a new Searcher object
    pub fn new(
        sa: Vec<i64>,
        sparseness_factor: u8,
        suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
        proteins: Proteins,
        taxon_id_calculator: TaxonAggregator,
        function_aggregator: FunctionAggregator
    ) -> Self {
        Self {
            sa,
            sparseness_factor,
            suffix_index_to_protein,
            proteins,
            taxon_id_calculator,
            function_aggregator
        }
    }
    
    /// Compares the `search_string` to the `suffix`
    /// During search this function performs extra logic since the suffix array is build with I == L, while ` self.proteins.input_string` is the original text where I != L
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide being searched in the suffix array
    /// * `suffix` - The current suffix from the suffix array we are comparing with in the binary search
    /// * `skip` - How many characters we can skip in the comparison because we already know these match
    /// * `bound` - Indicates if we are searching for the min of max bound
    ///
    /// # Returns
    ///
    /// The first argument is true if `bound` == `Minimum` and `search_string` <= `suffix` or if `bound` == `Maximum` and `search_string` >= `suffix`
    /// The second argument indicates how far the `suffix` and `search_string` matched
    fn compare(
        &self,
        search_string: &[u8],
        suffix: i64,
        skip: usize,
        bound: BoundSearch,
    ) -> (bool, usize) {
        let mut index_in_suffix = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;

        // Depending on if we are searching for the min of max bound our condition is different
        let condition_check = match bound {
            Minimum => |a: u8, b: u8| a < b,
            Maximum => |a: u8, b: u8| a > b,
        };

        // match as long as possible
        while index_in_search_string < search_string.len()
            && index_in_suffix < self.proteins.input_string.len()
            && (search_string[index_in_search_string]
                == self.proteins.input_string[index_in_suffix]
                || (search_string[index_in_search_string] == b'L'
                    && self.proteins.input_string[index_in_suffix] == b'I')
                || (search_string[index_in_search_string] == b'I'
                    && self.proteins.input_string[index_in_suffix] == b'L'))
        {
            index_in_suffix += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty search string should not be found)
        if !search_string.is_empty() {
            if index_in_search_string == search_string.len() {
                is_cond_or_equal = true
            } else if index_in_suffix < self.proteins.input_string.len() {
                // in our index every L was replaced by a I, so we need to replace them if we want to search in the right direction
                let peptide_char = if search_string[index_in_search_string] == b'L' {
                    b'I'
                } else {
                    search_string[index_in_search_string]
                };

                let protein_char = if self.proteins.input_string[index_in_suffix] == b'L' {
                    b'I'
                } else {
                    self.proteins.input_string[index_in_suffix]
                };

                is_cond_or_equal = condition_check(peptide_char, protein_char);
            }
        }

        (is_cond_or_equal, index_in_search_string)
    }
    
    /// Searches for the minimum or maximum bound for a string in the suffix array
    ///
    /// # Arguments
    /// * `bound` - Indicates if we are searching the minimum or maximum bound
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// The first argument is true if a match was found
    /// The second argument indicates the index of the minimum or maximum bound for the match (depending on `bound`)
    fn binary_search_bound(&self, bound: BoundSearch, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, bound);

            found |= lcp_center == search_string.len();

            // update the left and right bound, depending on if we are searching the min or max bound
            if retval && bound == Minimum || !retval && bound == Maximum {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }

        // handle edge case to search at index 0
        if right == 1 && left == 0 {
            let (retval, lcp_center) =
                self.compare(search_string, self.sa[0], min(lcp_left, lcp_right), bound);

            found |= lcp_center == search_string.len();

            if bound == Minimum && retval {
                right = 0;
            }
        }

        match bound {
            Minimum => (found, right),
            Maximum => (found, left),
        }
    }

    /// Searches for the minimum and maximum bound for a string in the suffix array
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    ///
    /// # Returns
    ///
    /// Returns the minimum and maximum bound of all matches in the suffix array, or `NoMatches` if no matches were found
    pub fn search_bounds(&self, search_string: &[u8]) -> BoundSearchResult {
        let (found_min, min_bound) = self.binary_search_bound(Minimum, search_string);

        if !found_min {
            return BoundSearchResult::NoMatches;
        }

        let (_, max_bound) = self.binary_search_bound(Maximum, search_string);

        BoundSearchResult::SearchResult((min_bound, max_bound + 1))
    }

    /// Searches for the suffixes matching a search string
    /// During search I and L can be equated
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide we are searching in the suffix array
    /// * `max_matches` - The maximum amount of matches processed, if more matches are found we don't process them
    /// * `equalize_i_and_l` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns all the matching suffixes
    #[inline]
    pub fn search_matching_suffixes(
        &self,
        search_string: &[u8],
        max_matches: usize,
        equalize_i_and_l: bool,
    ) -> SearchAllSuffixesResult {
        let mut matching_suffixes: Vec<i64> = vec![];
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        let mut skip: usize = 0;
        while skip < self.sparseness_factor as usize {
            let mut il_locations_start = 0;
            while il_locations_start < il_locations.len() && il_locations[il_locations_start] < skip {
                il_locations_start += 1;
            }
            let il_locations_current_suffix = &il_locations[il_locations_start..];
            let current_search_string_prefix = &search_string[..skip];
            let current_search_string_suffix = &search_string[skip..];
            let search_bound_result = self.search_bounds(&search_string[skip..]);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if let BoundSearchResult::SearchResult((min_bound, max_bound)) = search_bound_result {
                // try all the partially matched suffixes and store the matching suffixes in an array (stop when our max number of matches is reached)
                let mut sa_index = min_bound;
                while sa_index < max_bound {
                    let suffix = self.sa[sa_index] as usize;
                    // filter away matches where I was wrongfully equalized to L, and check the unmatched prefix
                    // when I and L equalized, we only need to check the prefix, not the whole match, when the prefix is 0, we don't need to check at all
                    if suffix >= skip
                        && ((skip == 0
                            || Self::check_prefix(
                        current_search_string_prefix,
                                &self.proteins.input_string[suffix - skip..suffix],
                                equalize_i_and_l,
                            ))
                            && Self::check_suffix(
                                skip,
                                il_locations_current_suffix,
                                current_search_string_suffix,
                                &self.proteins.input_string
                                    [suffix..suffix + search_string.len() - skip],
                                equalize_i_and_l,
                            ))
                    {
                        matching_suffixes.push((suffix - skip) as i64);

                        // return if max number of matches is reached
                        if matching_suffixes.len() >= max_matches {
                            return SearchAllSuffixesResult::MaxMatches(matching_suffixes);
                        }
                    }
                    sa_index += 1;
                }
            }
            skip += 1;
        }

        if matching_suffixes.is_empty() {
            SearchAllSuffixesResult::NoMatches
        } else {
            SearchAllSuffixesResult::SearchResult(matching_suffixes)
        }
    }

    /// Returns true of the prefixes are the same
    /// if `equalize_i_and_l` is set to true, L and I are considered the same
    ///
    /// # Arguments
    /// * `search_string_prefix` - The unchecked prefix of the string/peptide that is searched
    /// * `index_prefix` - The unchecked prefix from the protein from the suffix array
    /// * `equalize_i_and_l` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string_prefix` and `index_prefix` are considered the same, otherwise false
    #[inline]
    fn check_prefix(
        search_string_prefix: &[u8],
        index_prefix: &[u8],
        equalize_i_and_l: bool,
    ) -> bool {
        if equalize_i_and_l {
            search_string_prefix.iter().zip(index_prefix).all(
                |(&search_character, &index_character)| {
                    search_character == index_character
                        || (search_character == b'I' && index_character == b'L')
                        || (search_character == b'L' && index_character == b'I')
                },
            )
        } else {
            search_string_prefix == index_prefix
        }
    }

    /// Returns true of the search_string and index_string are equal
    /// This is automatically true if `equalize_i_and_l` is set to true, since there matched during search where I = L
    /// If `equalize_i_and_l` is set to false, we need to check if the I and L locations have the same character
    ///
    /// # Arguments
    /// * `skip` - The used skip factor during the search iteration
    /// * `il_locations` - The locations of the I's and L's in the **original** peptide
    /// * `search_string` - The peptide that is being searched, but already with the skipped prefix removed from it
    /// * `index_string` - The suffix that search_string matches with when I and L were equalized during search
    /// * `equalize_i_and_l` - True if we want to equate I and L during search, otherwise false
    ///
    /// # Returns
    ///
    /// Returns true if `search_string` and `index_string` are considered the same, otherwise false
    fn check_suffix(
        skip: usize,
        il_locations: &[usize],
        search_string: &[u8],
        index_string: &[u8],
        equalize_i_and_l: bool,
    ) -> bool {
        if equalize_i_and_l {
            true
        } else {
            for &il_location in il_locations {
                let index = il_location - skip;
                if search_string[index] != index_string[index] {
                    return false;
                }
            }
            true
        }
    }

    /// Returns all the proteins that correspond with the provided suffixes
    ///
    /// # Arguments
    /// * `suffixes` - List of suffix indices
    ///
    /// # Returns
    ///
    /// Returns the proteins that every suffix is a part of 
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &Vec<i64>) -> Vec<&Protein> {
        let mut res = vec![];
        for &suffix in suffixes {
            let protein_index = self.suffix_index_to_protein.suffix_to_protein(suffix);
            if !protein_index.is_null() {
                res.push(&self.proteins[protein_index as usize]);
            }
        }
        res
    }

    /// Searches all the matching proteins for a search_string/peptide in the suffix array
    ///
    /// # Arguments
    /// * `search_string` - The string/peptide being searched
    /// * `equalize_i_and_l` - If set to true, I and L are equalized during search
    ///
    /// # Returns
    ///
    /// Returns the matching proteins for the search_string
    pub fn search_proteins_for_peptide(&self, search_string: &[u8], equalize_i_and_l: bool) -> Vec<&Protein> {
        let mut matching_suffixes = vec![];
        if let SearchAllSuffixesResult::SearchResult(suffixes) =
            self.search_matching_suffixes(search_string, usize::MAX, equalize_i_and_l)
        {
            matching_suffixes = suffixes;
        }
        self.retrieve_proteins(&matching_suffixes)
    }

    /// Retrieves the taxonomic analysis for a collection of proteins
    ///
    /// # Arguments
    /// * `proteins` - A collection of proteins
    ///
    /// # Returns
    ///
    /// Returns the taxonomic analysis result for the given list of proteins
    #[inline]
    pub fn retrieve_lca(&self, proteins: &[&Protein]) -> Option<TaxonId> {
        let taxon_ids: Vec<TaxonId> = proteins.iter().map(|prot| prot.taxon_id).collect();
        
        self.taxon_id_calculator
            .aggregate(taxon_ids)
            .map(|id| self.taxon_id_calculator
                .snap_taxon(id)
            )
    }

    /// Returns true if the protein is considered valid by the provided taxonomy file
    ///
    /// # Arguments
    /// * `protein` - A protein of which we want to check the validity
    ///
    /// # Returns
    ///
    ///  Returns true if the protein is considered valid by the provided taxonomy file
    pub fn taxon_valid(&self, protein: &Protein) -> bool {
        self.taxon_id_calculator.taxon_valid(protein.taxon_id)
    }

    /// Retrieves the functional analysis for a collection of proteins
    ///
    /// # Arguments
    /// * `proteins` - A collection of proteins
    ///
    /// # Returns
    ///
    /// Returns the functional analysis result for the given list of proteins
    pub fn retrieve_function(&self, proteins: &[&Protein]) -> Option<FunctionalAggregation> {
        let res = self.function_aggregator.aggregate(proteins.to_vec());
        Some(res)
    }

    /// Retrieves the all the functional annotations for a collection of proteins
    ///
    /// # Arguments
    /// * `proteins` - A collection of proteins
    ///
    /// # Returns
    ///
    /// Returns all functional annotations for a collection of proteins
    pub fn get_all_functional_annotations(&self, proteins: &[&Protein]) -> Vec<Vec<String>> {
        self.function_aggregator.get_all_functional_annotations(proteins)
    }
    
}
