use crate::suffix_to_protein_index::SuffixToProteinIndex;
use crate::Nullable;
use std::cmp::min;
use std::collections::VecDeque;
use tsv_utils::taxon_id_calculator::TaxonIdCalculator;
use tsv_utils::{Protein, Proteins};
use umgap::taxon::TaxonId;

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

    fn compare(
        &self,
        search_string: &[u8],
        suffix: i64,
        skip: usize,
        compare_fn: fn(u8, u8) -> bool,
        equalize_i_and_l: bool,
    ) -> (bool, usize, Vec<usize>) {
        let mut il_locations = vec![];
        let mut index = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;
        // match as long as possible
        while index_in_search_string < search_string.len()
            && index < self.proteins.input_string.len()
            && search_string[index_in_search_string] == self.proteins.input_string[index]
        {
            // if we want to set I and L equal, we need to know where to "split"
            if equalize_i_and_l && (search_string[index_in_search_string] == b'I') {
                il_locations.push(index_in_search_string);
            }
            index += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically (and the empty search string should not be found)
        if !search_string.is_empty()
            && (index_in_search_string == search_string.len()
                || (index < self.proteins.input_string.len()
                    && compare_fn(
                        search_string[index_in_search_string],
                        self.proteins.input_string[index],
                    )))
        {
            is_cond_or_equal = true;
        }
        (is_cond_or_equal, index_in_search_string, il_locations)
    }

    fn binary_search_min_match(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
    ) -> (Vec<bool>, Vec<usize>) {
        let mut results = vec![];
        let mut found_array = vec![];
        let mut configurations_to_visit: VecDeque<(usize, usize, usize, usize, Vec<u8>)> =
            VecDeque::from([(0, self.sa.len(), 0, 0, search_string.to_owned())]);

        while !configurations_to_visit.is_empty() {
            let (mut left, mut right, mut lcp_left, mut lcp_right, search_string) =
                configurations_to_visit.pop_front().unwrap();
            let mut found = false;

            // repeat until search window is minimum size OR we matched the whole search string last iteration
            while right - left > 1 {
                let center = (left + right) / 2;
                let skip = min(lcp_left, lcp_right);
                let (retval, lcp_center, il_locations) = self.compare(
                    &search_string,
                    self.sa[center],
                    skip,
                    |a, b| a < b,
                    equalize_i_and_l,
                );
                if equalize_i_and_l {
                    for switch_location in il_locations {
                        let mut search_string_copy = search_string.clone();
                        search_string_copy[switch_location] =
                            match search_string_copy[switch_location] {
                                b'I' => b'L',
                                _ => search_string_copy[switch_location],
                            };

                        configurations_to_visit.push_back((
                            left,
                            right,
                            lcp_left,
                            lcp_right,
                            search_string_copy,
                        ))
                    }
                }

                found |= lcp_center == search_string.len();
                if retval {
                    right = center;
                    lcp_right = lcp_center;
                } else {
                    left = center;
                    lcp_left = lcp_center;
                }
            }

            // handle edge case to search at index 0
            if right == 1 && left == 0 {
                let (retval, lcp_center, il_locations) = self.compare(
                    &search_string,
                    self.sa[0],
                    min(lcp_left, lcp_right),
                    |a, b| a < b,
                    equalize_i_and_l,
                );
                if equalize_i_and_l {
                    for switch_location in il_locations {
                        let mut search_string_copy = search_string.to_owned();
                        search_string_copy[switch_location] =
                            match search_string_copy[switch_location] {
                                b'I' => b'L',
                                _ => search_string_copy[switch_location],
                            };

                        configurations_to_visit.push_back((
                            left,
                            right,
                            lcp_left,
                            lcp_right,
                            search_string_copy,
                        ))
                    }
                }

                found |= lcp_center == search_string.len();

                if retval {
                    right = 0;
                }
            }

            results.push(right);
            found_array.push(found)
        }

        (found_array, results)
    }

    fn binary_search_max_match(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
    ) -> (Vec<bool>, Vec<usize>) {
        let mut results = vec![];
        let mut found_array = vec![];
        let mut configurations_to_visit: VecDeque<(usize, usize, usize, usize, Vec<u8>)> =
            VecDeque::from([(0, self.sa.len(), 0, 0, search_string.to_owned())]);

        while !configurations_to_visit.is_empty() {
            let (mut left, mut right, mut lcp_left, mut lcp_right, search_string) =
                configurations_to_visit.pop_front().unwrap();
            let mut found = false;

            // repeat until search window is minimum size OR we matched the whole search string last iteration
            while right - left > 1 {
                let center = (left + right) / 2;
                let skip = min(lcp_left, lcp_right);
                let (retval, lcp_center, il_locations) = self.compare(
                    &search_string,
                    self.sa[center],
                    skip,
                    |a, b| a > b,
                    equalize_i_and_l,
                );

                if equalize_i_and_l {
                    for switch_location in il_locations {
                        let mut search_string_copy = search_string.to_owned();
                        search_string_copy[switch_location] =
                            match search_string_copy[switch_location] {
                                b'I' => b'L',
                                b'L' => b'I',
                                _ => search_string_copy[switch_location],
                            };

                        configurations_to_visit.push_back((
                            left,
                            right,
                            lcp_left,
                            lcp_right,
                            search_string_copy,
                        ))
                    }
                }

                found |= lcp_center == search_string.len();
                if retval {
                    left = center;
                    lcp_left = lcp_center;
                } else {
                    right = center;
                    lcp_right = lcp_center;
                }
            }

            results.push(left);
            found_array.push(found)
        }

        (found_array, results)
    }

    fn binary_search_match(&self, search_string: &[u8], equalize_i_and_l: bool) -> bool {
        let mut left: usize = 0;
        let mut right: usize = self.sa.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while !found && right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center, il_locations) = self.compare(
                search_string,
                self.sa[center],
                skip,
                |a, b| a < b,
                equalize_i_and_l,
            ); // TODO: update for equal I and L
            found = lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }

        // handle edge case to search at index 0
        if !found && right == 1 && left == 0 {
            let (_, lcp_center, il_locations) = self.compare(
                search_string,
                self.sa[0],
                min(lcp_left, lcp_right),
                |a, b| a < b,
                equalize_i_and_l,
            ); // TODO: update for equal I and L
            found = lcp_center == search_string.len();
        }

        found
    }

    pub fn search_bounds(
        &self,
        search_string: &[u8],
        equalize_i_and_l: bool,
    ) -> (bool, Vec<(usize, usize)>) {
        let (found, min_bound) = self.binary_search_min_match(search_string, equalize_i_and_l);
        if !found.iter().any(|&f| f) {
            return (false, vec![(0, self.sa.len())]);
        }
        let (_, max_bound) = self.binary_search_max_match(search_string, equalize_i_and_l);

        (
            true,
            found // zip the min and max bounds together, but filter out those who correspond to a not_found entry
                .iter()
                .zip(min_bound)
                .into_iter()
                .zip(max_bound.into_iter().map(|v| v + 1))
                .filter(|((&found, _), _)| found)
                .map(|((_, min_bound), max_bound)| (min_bound, max_bound))
                .collect(),
        ) // add 1 to max bound to have exclusive end
    }

    pub fn search_if_match(&self, search_string: &[u8], equalize_i_and_l: bool) -> bool {
        // case where SA is not sparse
        if self.sample_rate == 1 {
            return self.binary_search_match(search_string, equalize_i_and_l);
        }

        for skip in 0..self.sample_rate as usize {
            let (found, min_max_bounds) =
                self.search_bounds(&search_string[skip..], equalize_i_and_l);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if found {
                for (min_bound, max_bound) in min_max_bounds {
                    let unmatched_prefix = &search_string[..skip];
                    // try all the partially matched suffixes
                    for sa_index in min_bound..max_bound {
                        let suffix = self.sa[sa_index] as usize;
                        if suffix >= skip
                            && unmatched_prefix
                                == &self.proteins.input_string[suffix - skip..suffix]
                        {
                            return true;
                        }
                    }
                }
            }
        }
        false
    }

    /// Search all the suffixes that search string matches with
    /// The first value is a boolean indicating if the cutoff is used, the second value returns the actual taxa
    #[inline]
    pub fn search_matching_suffixes(
        &self,
        search_string: &[u8],
        max_matches: usize,
        equalize_i_and_l: bool,
    ) -> (bool, Vec<i64>) {
        let mut matching_suffixes: Vec<i64> = vec![];
        let mut skip: usize = 0;
        while skip < self.sample_rate as usize {
            let (found, min_max_bounds) =
                self.search_bounds(&search_string[skip..], equalize_i_and_l);
            // if the shorter part is matched, see if what goes before the matched suffix matches the unmatched part of the prefix
            if found {
                let unmatched_prefix = &search_string[..skip];
                for (min_bound, max_bound) in min_max_bounds {
                    // try all the partially matched suffixes and store the matching suffixes in an array (stop when our max number of matches is reached)
                    let mut sa_index = min_bound;
                    while sa_index < max_bound {
                        let suffix = self.sa[sa_index] as usize;
                        // if skip is 0, then we already checked the complete match during bound search, otherwise check if the skipped part also matches
                        if skip == 0
                            || (suffix >= skip
                                && unmatched_prefix
                                    == &self.proteins.input_string[suffix - skip..suffix])
                        {
                            matching_suffixes.push((suffix - skip) as i64);

                            // return if max number of matches is reached
                            if matching_suffixes.len() >= max_matches {
                                return (true, matching_suffixes);
                            }
                        }

                        sa_index += 1;
                    }
                }
            }
            skip += 1;
        }
        (false, matching_suffixes)
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
        let (_, matching_suffixes) =
            self.search_matching_suffixes(search_string, usize::MAX, equalize_i_and_l);
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
