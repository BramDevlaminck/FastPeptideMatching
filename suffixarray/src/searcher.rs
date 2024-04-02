use std::cmp::min;
use std::collections::VecDeque;
use std::fmt::{Display, Formatter};

use umgap::taxon::TaxonId;

use tsv_utils::taxon_id_calculator::TaxonIdCalculator;
use tsv_utils::{Protein, Proteins};

use crate::searcher::BoundSearch::{Maximum, Minimum};
use crate::searcher::SearchMatchResult::{Found, NotFound, OutOfTime};
use crate::sequence_bitpattern::SequenceBitPattern;
use crate::suffix_to_protein_index::SuffixToProteinIndex;
use crate::util::get_time_ms;
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
    OutOfTime,
    NoMatches,
    SearchResult(Vec<(usize, usize)>),
}

#[derive(Debug)]
pub enum SearchAllSuffixesResult {
    OutOfTime,
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
            (SearchAllSuffixesResult::OutOfTime, SearchAllSuffixesResult::OutOfTime) => true,
            (
                SearchAllSuffixesResult::MaxMatches(arr1),
                SearchAllSuffixesResult::MaxMatches(arr2),
            ) => array_eq_unordered(arr1, arr2),
            (
                SearchAllSuffixesResult::SearchResult(arr1),
                SearchAllSuffixesResult::SearchResult(arr2),
            ) => array_eq_unordered(arr1, arr2),
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
            Found => write!(f, "true")
        }
    }
}

pub struct Searcher {
    sa: Vec<i64>,
    pub sample_rate: u8,
    suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
    proteins: Proteins,
    taxon_id_calculator: TaxonIdCalculator,
}

impl Searcher {
    pub fn new(
        sa: Vec<i64>,
        sample_rate: u8,
        suffix_index_to_protein: Box<dyn SuffixToProteinIndex>,
        proteins: Proteins,
        taxon_id_calculator: TaxonIdCalculator,
    ) -> Self {
        Self {
            sa,
            sample_rate,
            suffix_index_to_protein,
            proteins,
            taxon_id_calculator,
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
        equalize_i_and_l: bool,
        bound: BoundSearch,
    ) -> (bool, usize, Vec<usize>) {
        let mut positions_to_switch_i_and_l = vec![];
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
            && search_string[index_in_search_string] == self.proteins.input_string[index_in_suffix]
        {
            // if we want to set I and L equal, we need to know where to "split"
            if equalize_i_and_l && search_string[index_in_search_string] == b'I'
            {
                positions_to_switch_i_and_l.push(index_in_search_string);
            }
            index_in_suffix += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty search string should not be found)
        if !search_string.is_empty()
            && (index_in_search_string == search_string.len()
                || (index_in_suffix < self.proteins.input_string.len()
                    && condition_check(
                        search_string[index_in_search_string],
                        self.proteins.input_string[index_in_suffix],
                    )))
        {
            is_cond_or_equal = true;
        }

        // Pattern of searching an 'L' is only different from searching a 'I' in 3 cases. This is when the character we compare it with is a 'J', 'K' or 'L'
        if equalize_i_and_l
            && index_in_search_string < search_string.len()
            && index_in_suffix < self.proteins.input_string.len()
            && search_string[index_in_search_string] == b'I'
            && (self.proteins.input_string[index_in_suffix] > b'I'
                && self.proteins.input_string[index_in_suffix] <= b'L')
        {
            positions_to_switch_i_and_l.push(index_in_search_string);
        }

        (
            is_cond_or_equal,
            index_in_search_string,
            positions_to_switch_i_and_l,
        )
    }

    /// When `equalize_i_and_l` is set, this function adds the needed variants to the `to_visit` deque.
    /// This is done using the correct bounds (`left`, `right`, `lcp_left` and `lcp_right`).
    /// To prevent searching the same string multiple times we use a kind of bitvector (the `visited_strings` variable to keep track of which variants of the input sequence are already in the queue to be searched)
    /// The `il_locations` vector contains all indices where there is an I or L in the `current_search_string`
    fn replace_i_with_l(
        equalize_i_and_l: bool,
        to_visit: &mut VecDeque<(usize, usize, usize, usize, Vec<u8>)>,
        current_search_string: &mut [u8],
        positions_to_switch_i_and_l: Vec<usize>,
        left: usize,
        right: usize,
        lcp_left: usize,
        lcp_right: usize,
        visited_strings: &mut SequenceBitPattern,
        il_locations: &Vec<usize>,
    ) {
        if equalize_i_and_l {
            for switch_location in positions_to_switch_i_and_l {
                let mut search_string_copy = current_search_string.to_owned();
                search_string_copy[switch_location] = match search_string_copy[switch_location] {
                    b'I' => b'L',
                    _ => search_string_copy[switch_location],
                };

                // only add the IL variant of this string if it is not yet visited (or in the queue waiting to be visited)
                if !visited_strings.check_if_contains_and_add(&search_string_copy, il_locations) {
                    to_visit.push_back((
                        left,
                        right,
                        lcp_left,
                        lcp_right,
                        search_string_copy.clone(),
                    ));
                }
            }
        }
    }

    fn binary_search_bound(
        &self,
        bound: BoundSearch,
        search_string: &[u8],
        il_locations: &Vec<usize>,
        end_search_time: f64,
    ) -> BinarySearchBoundResult {
        let equalize_i_and_l = !il_locations.is_empty();
        let mut visited_pattern = SequenceBitPattern::new(il_locations);
        let mut results = vec![];
        let mut found_array = vec![];
        let mut configurations_to_visit: VecDeque<(usize, usize, usize, usize, Vec<u8>)> =
            VecDeque::from([(0, self.sa.len(), 0, 0, search_string.to_owned())]);

        while let Some((mut left, mut right, mut lcp_left, mut lcp_right, mut search_string)) =
            configurations_to_visit.pop_front()
        {
            if get_time_ms().unwrap() > end_search_time {
                return BinarySearchBoundResult::OutOfTime
            }

            let mut found = false;

            // repeat until search window is minimum size OR we matched the whole search string last iteration
            while right - left > 1 {
                let center = (left + right) / 2;
                let skip = min(lcp_left, lcp_right);
                let (retval, lcp_center, positions_to_switch_i_and_l) = self.compare(
                    &search_string,
                    self.sa[center],
                    skip,
                    equalize_i_and_l,
                    bound,
                );

                Self::replace_i_with_l(
                    equalize_i_and_l,
                    &mut configurations_to_visit,
                    &mut search_string,
                    positions_to_switch_i_and_l,
                    left,
                    right,
                    lcp_left,
                    lcp_right,
                    &mut visited_pattern,
                    il_locations,
                );

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
                let (retval, lcp_center, positions_to_switch_i_and_l) = self.compare(
                    &search_string,
                    self.sa[0],
                    min(lcp_left, lcp_right),
                    equalize_i_and_l,
                    bound,
                );

                // handle I's we passed while progressing
                Self::replace_i_with_l(
                    equalize_i_and_l,
                    &mut configurations_to_visit,
                    &mut search_string,
                    positions_to_switch_i_and_l,
                    left,
                    right,
                    lcp_left,
                    lcp_right,
                    &mut visited_pattern,
                    il_locations,
                );

                found |= lcp_center == search_string.len();

                if bound == Minimum && retval {
                    right = 0;
                }
            }

            match bound {
                Minimum => results.push(right),
                Maximum => results.push(left),
            }
            found_array.push(found)
        }

        BinarySearchBoundResult::SearchResult((found_array, results))
    }

    /// Search the min and max bound of the `search_string` in the suffix array.
    /// If `equalize_i_and_l` is set to true we will also allow an L on every location an I is present.
    /// This means that we expect the `search_string` to be preprocessed when we want to equalize I and L.
    /// Every L should already be replaced by an I in the `search_string`!
    pub fn search_bounds(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
        end_search_time: f64,
    ) -> BoundSearchResult {
        // if we equalize I and L, calculate the locations of these letters, otherwise just use the empty vector, since it will be unused
        // we will use the size as this vector to indicate if we should equalize I and L or not.
        // it is possible that equalize_i_and_l was set to true, but no I or L is part of the string, in that case we will spare all the extra checks if were are at an I or L
        let mut il_locations = vec![];
        if equalize_i_and_l {
            search_string
                .iter()
                .enumerate()
                .for_each(|(index, &character)| {
                    if character == b'I' || character == b'L' {
                        il_locations.push(index)
                    }
                });
            // 2^34 options in a bitvector is around 2 GB. Going allowing even more I's and L's will allocate a huge amount of memory
            if il_locations.len() > 34 {
                return BoundSearchResult::OutOfTime;
            }
        }
        


        let (found_min, min_bound) = match self.binary_search_bound(
            Minimum,
            search_string,
            &il_locations,
            end_search_time,
        ) {
            BinarySearchBoundResult::OutOfTime => return BoundSearchResult::OutOfTime,
            BinarySearchBoundResult::SearchResult((found_min, found_max)) => (found_min, found_max),
        };

        if !found_min.iter().any(|&f| f) {
            return BoundSearchResult::NoMatches;
        }

        let (found_max, max_bound) = match self.binary_search_bound(
            Maximum,
            search_string,
            &il_locations,
            end_search_time,
        ) {
            BinarySearchBoundResult::OutOfTime => return BoundSearchResult::OutOfTime,
            BinarySearchBoundResult::SearchResult((found_min, found_max)) => (found_min, found_max),
        };

        // Only get the values from the min and max bound search that actually had a match
        let min_bounds: Vec<usize> = found_min
            .iter()
            .zip(min_bound)
            .filter(|(&found, _)| found)
            .map(|(_, bound)| bound)
            .collect();
        let max_bounds: Vec<usize> = found_max
            .iter()
            .zip(max_bound)
            .filter(|(&found, _)| found)
            .map(|(_, bound)| bound + 1)
            .collect();

        BoundSearchResult::SearchResult(min_bounds.into_iter().zip(max_bounds).collect())
    }

    /// Search if there is a match for the `search_string` can be found somewhere in the suffix array.
    /// If `equalize_i_and_l` is set to true we will also allow an L on every location an I is present.
    /// This means that we expect the `search_string` to be preprocessed when we want to equalize I and L.
    /// Every L should already be replaced by an I in the `search_string`!
    pub fn search_if_match(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
        max_peptide_search_time: MaxPeptideSearchTime,
    ) -> SearchMatchResult {
        let end_search_time: f64 =
            get_time_ms().unwrap() + Into::<f64>::into(max_peptide_search_time);

        for skip in 0..self.sample_rate as usize {
            let bound_search_res =
                self.search_bounds(&search_string[skip..], equalize_i_and_l, end_search_time);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if let BoundSearchResult::SearchResult(min_max_bounds) = bound_search_res {
                for (min_bound, max_bound) in min_max_bounds {
                    let unmatched_prefix = &search_string[..skip];
                    // try all the partially matched suffixes
                    for sa_index in min_bound..max_bound {
                        let suffix = self.sa[sa_index] as usize;
                        if suffix >= skip
                            && Self::equals_using_il_equality(
                                unmatched_prefix,
                                &self.proteins.input_string[suffix - skip..suffix],
                                equalize_i_and_l,
                            )
                        {
                            return Found;
                        }
                    }
                }
            } else if let BoundSearchResult::OutOfTime = bound_search_res {
                return OutOfTime;
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
        max_peptide_search_time: MaxPeptideSearchTime,
    ) -> SearchAllSuffixesResult {
        let end_search_time: f64 =
            get_time_ms().unwrap() + Into::<f64>::into(max_peptide_search_time);

        let mut matching_suffixes: Vec<i64> = vec![];
        let mut skip: usize = 0;
        while skip < self.sample_rate as usize {
            let search_bound_result =
                self.search_bounds(&search_string[skip..], equalize_i_and_l, end_search_time);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if let BoundSearchResult::SearchResult(min_max_bounds) = search_bound_result {
                let unmatched_prefix = &search_string[..skip];
                for (min_bound, max_bound) in min_max_bounds {
                    // try all the partially matched suffixes and store the matching suffixes in an array (stop when our max number of matches is reached)
                    let mut sa_index = min_bound;
                    while sa_index < max_bound {
                        let suffix = self.sa[sa_index] as usize;
                        // if skip is 0, then we already checked the complete match during bound search, otherwise check if the skipped part also matches
                        if skip == 0
                            || (suffix >= skip
                                && Self::equals_using_il_equality(
                                    unmatched_prefix,
                                    &self.proteins.input_string[suffix - skip..suffix],
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
            } else if let BoundSearchResult::OutOfTime = search_bound_result {
                return SearchAllSuffixesResult::OutOfTime;
            }
            skip += 1;
        }
        SearchAllSuffixesResult::SearchResult(matching_suffixes)
    }

    /// Returns true of the prefixes are the same
    /// if `equalize_i_and_l` is set to true, L and I are considered the same
    /// This implementation abuses the fact that every L is translated to a I in the search_prefix when we want to do this
    #[inline]
    fn equals_using_il_equality(
        search_string_prefix: &[u8],
        index_prefix: &[u8],
        equalize_i_and_l: bool,
    ) -> bool {
        if equalize_i_and_l {
            search_string_prefix.iter().zip(index_prefix).all(
                |(&search_character, &index_character)| {
                    search_character == index_character
                        || (search_character == b'I' && index_character == b'L')
                },
            )
        } else {
            search_string_prefix == index_prefix
        }
    }

    /// get all the proteins matching with the given suffixes
    #[inline]
    pub fn retrieve_proteins(&self, suffixes: &Vec<i64>) -> Vec<&Protein> {
        let mut res = vec![];
        for &suffix in suffixes {
            let protein_index = self.suffix_index_to_protein.suffix_to_protein(suffix);
            if !protein_index.is_null() {
                res.push(&self.proteins.proteins[protein_index as usize]);
            }
        }
        res
    }

    /// Search all the Proteins that a given search_string matches with
    pub fn search_protein(&self, search_string: &[u8], equalize_i_and_l: bool) -> Vec<&Protein> {
        let mut matching_suffixes = vec![];
        if let SearchAllSuffixesResult::SearchResult(suffixes) = self.search_matching_suffixes(
            search_string,
            usize::MAX,
            equalize_i_and_l,
            MaxPeptideSearchTime::Unlimited,
        ) {
            matching_suffixes = suffixes;
        }
        self.retrieve_proteins(&matching_suffixes)
    }

    #[inline]
    pub fn retrieve_lca(&self, proteins: &[&Protein]) -> Option<TaxonId> {
        let taxon_ids: Vec<TaxonId> = proteins.iter().map(|prot| prot.id).collect();
        match taxon_ids.is_empty() {
            true => None,
            false => Some(
                self.taxon_id_calculator
                    .snap_taxon_id(self.taxon_id_calculator.get_aggregate(taxon_ids)),
            ),
        }
    }

    /// Fetch the UniProt accession and taxa for all the proteins
    pub fn get_uniprot_and_taxa_ids(proteins: &[&Protein]) -> (Vec<String>, Vec<TaxonId>) {
        let mut uniprot_ids = vec![];
        let mut taxa = vec![];

        for &protein in proteins {
            taxa.push(protein.id);
            uniprot_ids.push(protein.uniprot_id.clone());
        }

        (uniprot_ids, taxa)
    }

    pub fn search_lca(&self, search_string: &[u8], equalize_i_and_l: bool) -> Option<TaxonId> {
        self.retrieve_lca(&self.search_protein(search_string, equalize_i_and_l))
    }
}

#[cfg(test)]
mod tests {
    use tsv_utils::taxon_id_calculator::{AggregationMethod, TaxonIdCalculator};
    use tsv_utils::{Protein, Proteins};

    use crate::searcher::{
        BoundSearchResult, MaxPeptideSearchTime, SearchAllSuffixesResult, SearchMatchResult,
        Searcher,
    };
    use crate::suffix_to_protein_index::SparseSuffixToProtein;

    fn get_example_proteins() -> Proteins {
        let text = "AI-BLACVAA-AC-KCRLZ$".to_string().into_bytes();
        Proteins {
            input_string: text,
            proteins: vec![
                Protein {
                    uniprot_id: String::new(),
                    sequence: (0, 2),
                    id: 0,
                },
                Protein {
                    uniprot_id: String::new(),
                    sequence: (3, 10),
                    id: 0,
                },
                Protein {
                    uniprot_id: String::new(),
                    sequence: (11, 13),
                    id: 0,
                },
                Protein {
                    uniprot_id: String::new(),
                    sequence: (14, 29),
                    id: 0,
                },
            ],
        }
    }

    #[test]
    fn test_search_simple() {
        let proteins = get_example_proteins();
        let sa = vec![
            19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 14, 4, 17, 16, 7, 18,
        ];

        let searcher = Searcher::new(
            sa,
            1,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        let end_time = f64::INFINITY;

        // search bounds 'A'
        let bounds_res = searcher.search_bounds(&[b'A'], false, end_time);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult(vec![(4, 9)]));

        // search bounds '$'
        let bounds_res = searcher.search_bounds(&[b'$'], false, end_time);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult(vec![(0, 1)]));

        // search bounds 'AC'
        let bounds_res = searcher.search_bounds(&[b'A', b'C'], false, end_time);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult(vec![(6, 8)]));
    }

    #[test]
    fn test_search_sparse() {
        let proteins = get_example_proteins();
        let sa = vec![9, 0, 3, 12, 15, 6, 18];

        let searcher = Searcher::new(
            sa,
            3,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search suffix 'VAA'
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'V', b'A', b'A'],
            usize::MAX,
            false,
            MaxPeptideSearchTime::Default,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![7])
        );

        // search suffix 'AC'
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'A', b'C'],
            usize::MAX,
            false,
            MaxPeptideSearchTime::Default,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![5, 11])
        );
    }

    #[test]
    fn test_il_equality() {
        let proteins = get_example_proteins();
        let sa = vec![
            19, 10, 2, 13, 9, 8, 11, 5, 0, 3, 12, 15, 6, 1, 14, 4, 17, 16, 7, 18,
        ];

        let searcher = Searcher::new(
            sa,
            1,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        let end_time = f64::INFINITY;

        // search bounds 'I' with equal I and L
        let bounds_res = searcher.search_bounds(&[b'I'], true, end_time);
        assert_eq!(
            bounds_res,
            BoundSearchResult::SearchResult(vec![(13, 14), (15, 17)])
        );

        // search bounds 'RIZ' with equal I and L
        let bounds_res = searcher.search_bounds(&[b'R', b'I', b'Z'], true, end_time);
        assert_eq!(bounds_res, BoundSearchResult::SearchResult(vec![(17, 18)]));
    }

    #[test]
    fn test_il_equality_sparse() {
        let proteins = get_example_proteins();
        let sa = vec![9, 0, 3, 12, 15, 6, 18];

        let searcher = Searcher::new(
            sa,
            3,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search bounds 'RIZ' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'R', b'I', b'Z'],
            usize::MAX,
            true,
            MaxPeptideSearchTime::Unlimited,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![16])
        );
    }

    // test edge case where an I or L is the first index in the sparse SA.
    #[test]
    fn test_l_first_index_in_sa() {
        let text = "LMOXZ$".to_string().into_bytes();

        let proteins = Proteins {
            input_string: text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                sequence: (0, 6),
                id: 0,
            }],
        };

        let sparse_sa = vec![0, 2, 4];
        let searcher = Searcher::new(
            sparse_sa,
            2,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'I', b'M'],
            usize::MAX,
            true,
            MaxPeptideSearchTime::Unlimited,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![0])
        );
    }

    #[test]
    fn test_matches() {
        let proteins = get_example_proteins();
        let sa = vec![10, 2, 8, 0, 12, 6, 14, 4, 16, 18];

        let searcher = Searcher::new(
            sa,
            2,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        assert_eq!(
            searcher.search_if_match(&[b'A', b'I'], false, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::Found
        );
        assert_eq!(
            searcher.search_if_match(&[b'B', b'L'], false, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::Found
        );
        assert_eq!(
            searcher.search_if_match(&[b'K', b'C', b'R'], false, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::Found
        );
        assert_eq!(
            searcher.search_if_match(&[b'A', b'L'], false, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::NotFound
        );
        assert_eq!(
            searcher.search_if_match(&[b'B', b'I'], false, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::NotFound
        );
        assert_eq!(
            searcher.search_if_match(&[b'A', b'I'], true, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::Found
        );
        assert_eq!(
            searcher.search_if_match(&[b'B', b'I'], true, MaxPeptideSearchTime::Unlimited),
            SearchMatchResult::Found
        );
        assert_eq!(
            searcher.search_if_match(
                &[b'K', b'C', b'R', b'I'],
                true,
                MaxPeptideSearchTime::Unlimited
            ),
            SearchMatchResult::Found
        );
    }

    #[test]
    fn test_il_missing_matches() {
        let text = "AAILLL$".to_string().into_bytes();

        let proteins = Proteins {
            input_string: text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                sequence: (0, 13),
                id: 0,
            }],
        };

        let sparse_sa = vec![6, 0, 1, 2, 5, 4, 3];
        let searcher = Searcher::new(
            sparse_sa,
            1,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'I'],
            usize::MAX,
            true,
            MaxPeptideSearchTime::Unlimited,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![2, 3, 4, 5])
        );
    }

    #[test]
    fn test_il_duplication() {
        let text = "IIIILL$".to_string().into_bytes();

        let proteins = Proteins {
            input_string: text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                sequence: (0, 13),
                id: 0,
            }],
        };

        let sparse_sa = vec![6, 0, 1, 2, 3, 5, 4];
        let searcher = Searcher::new(
            sparse_sa,
            1,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'I', b'I'],
            usize::MAX,
            true,
            MaxPeptideSearchTime::Unlimited,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4])
        );
    }

    #[test]
    fn test_il_duplication2() {
        let text = "IILLLL$".to_string().into_bytes();

        let proteins = Proteins {
            input_string: text,
            proteins: vec![Protein {
                uniprot_id: String::new(),
                sequence: (0, 13),
                id: 0,
            }],
        };

        let sparse_sa = vec![6, 0, 1, 5, 4, 3, 2];
        let searcher = Searcher::new(
            sparse_sa,
            1,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search bounds 'IM' with equal I and L
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'I', b'I'],
            usize::MAX,
            true,
            MaxPeptideSearchTime::Unlimited,
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::SearchResult(vec![0, 1, 2, 3, 4])
        );
    }

    #[test]
    fn test_timeout() {
        let proteins = get_example_proteins();
        let sa = vec![9, 0, 3, 12, 15, 6, 18];

        let searcher = Searcher::new(
            sa,
            3,
            Box::new(SparseSuffixToProtein::new(&proteins.input_string)),
            proteins,
            *TaxonIdCalculator::new("../small_taxonomy.tsv", AggregationMethod::LcaStar),
        );

        // search suffix 'VAA'
        let found_suffixes = searcher.search_matching_suffixes(
            &[b'V', b'A', b'A'],
            usize::MAX,
            false,
            MaxPeptideSearchTime::Value(0.0),
        );
        assert_eq!(
            found_suffixes,
            SearchAllSuffixesResult::OutOfTime
        );

        // search suffix 'AC'
        let found_suffixes = searcher.search_if_match(
            &[b'A', b'C'],
            false,
            MaxPeptideSearchTime::Value(0.0),
        );
        assert_eq!(
            found_suffixes,
            SearchMatchResult::OutOfTime
        );
    }
}
