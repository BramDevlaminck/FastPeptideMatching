use std::cmp::min;
use std::fmt::{Display, Formatter};

use bitarray::BitArray;
use sa_mappings::functionality::{FunctionAggregator, FunctionalAggregation};
use sa_mappings::proteins::{Protein, Proteins};
use sa_mappings::taxonomy::TaxonAggregator;
use umgap::taxon::TaxonId;

use crate::sa_searcher::BoundSearch::{Maximum, Minimum};
use crate::sa_searcher::SearchMatchResult::{Found, NotFound, OutOfTime};
use crate::suffix_to_protein_index::SuffixToProteinIndex;
use crate::Nullable;

/// Enum indicating if we are searching for the minimum, or maximum bound in the suffix array
#[derive(Clone, Copy, PartialEq)]
enum BoundSearch {
    Minimum,
    Maximum,
}

/// Enum to how many ms we want to spend at most to search 1 peptide
pub enum MaxPeptideSearchTime {
    Default,
    Value(f64),
    Unlimited,
}

impl Into<f64> for MaxPeptideSearchTime {
    fn into(self) -> f64 {
        match self {
            MaxPeptideSearchTime::Default => 60000.0,
            MaxPeptideSearchTime::Value(search_time) => search_time,
            MaxPeptideSearchTime::Unlimited => f64::INFINITY,
        }
    }
}

// Enums to present the possible results for different functions

#[derive(PartialEq, Debug)]
pub enum BinarySearchBoundResult {
    OutOfTime,
    SearchResult((Vec<bool>, Vec<usize>)),
}

#[derive(PartialEq, Debug)]
pub enum BoundSearchResult {
    NoMatches,
    SearchResult((usize, usize)),
}

#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    NoMatches,
    MaxMatches(Vec<i64>),
    SearchResult(Vec<i64>),
}

impl PartialEq for SearchAllSuffixesResult {
    fn eq(&self, other: &Self) -> bool {
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

#[derive(PartialEq, Debug)]
pub enum SearchMatchResult {
    OutOfTime,
    NotFound,
    Found,
}

impl Display for SearchMatchResult {
    fn fmt(&self, f: &mut Formatter<'_>) -> std::fmt::Result {
        match self {
            OutOfTime => write!(f, "search time limit reached"),
            NotFound => write!(f, "false"),
            Found => write!(f, "true"),
        }
    }
}

pub struct Searcher {
    sa: BitArray<38>,
    pub sample_rate: u8,
    suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
    proteins: Proteins,
    taxon_id_calculator: TaxonAggregator,
    function_aggregator: FunctionAggregator
}

impl Searcher {
    pub fn new(
        sa: BitArray<38>,
        sample_rate: u8,
        suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
        proteins: Proteins,
        taxon_id_calculator: TaxonAggregator,
        function_aggregator: FunctionAggregator
    ) -> Self {
        Self {
            sa,
            sample_rate,
            suffix_index_to_protein,
            proteins,
            taxon_id_calculator,
            function_aggregator
        }
    }

    /// Compare function to perform the binary search
    /// This function forwards progresses the search of the `search_string` as much as possible on the current `suffix`.
    /// This search is a bit different depending on the current bound (searching for the min or max)
    ///
    /// When `equalize_i_and_l` is set to true this function also keeps track of which I's and L's we have progressed through and if we should start a new branch in the search tree
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
            let (retval, lcp_center) = self.compare(search_string, self.sa.get(center) as i64, skip, bound);

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
                self.compare(search_string, self.sa.get(0) as i64, min(lcp_left, lcp_right), bound);

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

    /// Search the min and max bound of the `search_string` in the suffix array.
    /// If `equalize_i_and_l` is set to true we will also allow an L on every location an I is present.
    /// This means that we expect the `search_string` to be preprocessed when we want to equalize I and L.
    /// Every L should already be replaced by an I in the `search_string`!
    pub fn search_bounds(&self, search_string: &[u8]) -> BoundSearchResult {
        let (found_min, min_bound) = self.binary_search_bound(Minimum, search_string);

        if !found_min {
            return BoundSearchResult::NoMatches;
        }

        let (_, max_bound) = self.binary_search_bound(Maximum, search_string);

        BoundSearchResult::SearchResult((min_bound, max_bound + 1))
    }

    /// Search if there is a match for the `search_string` can be found somewhere in the suffix array.
    /// If `equalize_i_and_l` is set to true we will also allow an L on every location an I is present.
    /// This means that we expect the `search_string` to be preprocessed when we want to equalize I and L.
    /// Every L should already be replaced by an I in the `search_string`!
    pub fn search_if_match(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
    ) -> SearchMatchResult {
        let mut il_locations = vec![];
        for (i, &character) in search_string.iter().enumerate() {
            if character == b'I' || character == b'L' {
                il_locations.push(i);
            }
        }

        for skip in 0..self.sample_rate as usize {
            let mut il_locations_start = 0;
            while il_locations_start < il_locations.len() && il_locations[il_locations_start] < skip {
                il_locations_start += 1;
            }
            let il_locations_current_suffix = &il_locations[il_locations_start..];
            let current_search_string_prefix = &search_string[..skip];
            let current_search_string_suffix = &search_string[skip..];
            let bound_search_res = self.search_bounds(&search_string[skip..]);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if let BoundSearchResult::SearchResult((min_bound, max_bound)) = bound_search_res {
                // try all the partially matched suffixes
                for sa_index in min_bound..max_bound {
                    let suffix = self.sa.get(sa_index) as usize;
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
                        return Found;
                    }
                }
            }
        }
        NotFound
    }

    /// Search all the suffixes that search string matches with
    /// The first value is a boolean indicating if the cutoff is used, the second value returns the actual taxa
    ///
    /// If `equalize_i_and_l` is set to true we will also allow an L on every location an I is present.
    /// This means that we expect the `search_string` to be preprocessed when we want to equalize I and L.
    /// Every L should already be replaced by an I in the `search_string`!
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
        while skip < self.sample_rate as usize {
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
                    let suffix = self.sa.get(sa_index) as usize;
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

    /// Returns true if I or L equality is used, since the suffix will surely match then
    /// otherwise we have to check the locations where an I or L was present
    /// If there is a mismatch on 1 of these locations, the suffix does not match
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

    /// get all the proteins matching with the given suffixes
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

    /// Search all the Proteins that a given search_string matches with
    pub fn search_protein(&self, search_string: &[u8], equalize_i_and_l: bool) -> Vec<&Protein> {
        let mut matching_suffixes = vec![];
        if let SearchAllSuffixesResult::SearchResult(suffixes) =
            self.search_matching_suffixes(search_string, usize::MAX, equalize_i_and_l)
        {
            matching_suffixes = suffixes;
        }
        self.retrieve_proteins(&matching_suffixes)
    }

    #[inline]
    pub fn retrieve_lca(&self, proteins: &[&Protein]) -> Option<TaxonId> {
        let taxon_ids: Vec<TaxonId> = proteins.iter().map(|prot| prot.taxon_id).collect();
        
        self.taxon_id_calculator
            .aggregate(taxon_ids)
            .map(|id| self.taxon_id_calculator
                .snap_taxon(id)
            )
    }
    
    /// Return whether the taxon id of the protein has its valid bit set in the taxonomy
    pub fn taxon_valid(&self, protein: &Protein) -> bool {
        self.taxon_id_calculator.taxon_valid(protein.taxon_id)
    }
    
    pub fn retrieve_function(&self, proteins: &[&Protein]) -> Option<FunctionalAggregation> {
        let res = self.function_aggregator.aggregate(proteins.to_vec());
        Some(res)
    }

    /// Fetch the UniProt accession and taxa for all the proteins
    pub fn get_uniprot_and_taxa_ids(proteins: &[&Protein]) -> (Vec<String>, Vec<TaxonId>) {
        let mut uniprot_ids = vec![];
        let mut taxa = vec![];

        for &protein in proteins {
            taxa.push(protein.taxon_id);
            uniprot_ids.push(protein.uniprot_id.clone());
        }

        (uniprot_ids, taxa)
    }
    
}