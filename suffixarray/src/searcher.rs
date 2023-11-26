use std::cmp::{max, min};
use tsv_utils::Protein;

pub struct Searcher<'a> {
    original_input_string: &'a [u8],
    sa: &'a Vec<i32>,
    proteins: &'a Vec<Protein>,
}

impl <'a> Searcher<'a> {

    pub fn new(original_input_string: &'a [u8], sa: &'a Vec<i32>, proteins: &'a Vec<Protein>) -> Self {
        Self {
            original_input_string,
            sa,
            proteins
        }
    }

    fn compare(&self, search_string: &[u8], suffix: i32, skip: usize, compare_fn: fn(usize, usize) -> bool) -> (bool, usize) {
        let mut index = (suffix as usize) + skip;
        let mut index_in_search_string = skip;
        let mut is_cond_or_equal = false;
        // match as long as possible
        while index_in_search_string < search_string.len()
            && search_string[index_in_search_string] == self.original_input_string[index] {
            index += 1;
            index_in_search_string += 1;
        }
        // check if match found OR current search string is smaller lexicographically
        if index_in_search_string == search_string.len() || compare_fn(search_string[index_in_search_string] as usize, self.original_input_string[index] as usize) {
            is_cond_or_equal = true;
        }
        (is_cond_or_equal, index_in_search_string)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_min_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a < b);
            found |= lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }
        (found, right)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_max_match(&self, search_string: &[u8]) -> (bool, usize) {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a > b);
            found |= lcp_center == search_string.len();
            if retval {
                left = center;
                lcp_left = lcp_center;
            } else {
                right = center;
                lcp_right = lcp_center;
            }
        }
        (found, left)
    }

    // note: $ cannot be found since this is the first value in the array and the way the bounds are changed with center, index 0 is never checked
    fn binary_search_match(&self, search_string: &[u8]) -> bool {
        let mut left: usize = 0;
        let mut right: usize = self.original_input_string.len();
        let mut lcp_left: usize = 0;
        let mut lcp_right: usize = 0;
        let mut found = false;

        // repeat until search window is minimum size OR we matched the whole search string last iteration
        while !found && right - left > 1 {
            let center = (left + right) / 2;
            let skip = min(lcp_left, lcp_right);
            let (retval, lcp_center) = self.compare(search_string, self.sa[center], skip, |a, b| a < b);
            found |= lcp_center == search_string.len();
            if retval {
                right = center;
                lcp_right = lcp_center;
            } else {
                left = center;
                lcp_left = lcp_center;
            }
        }
        found
    }

    pub fn search_bounds(&self, search_string: &[u8]) -> (bool, usize, usize) {
        let (found, min_bound) = self.binary_search_min_match(search_string);
        if !found {
            return (false, 0, self.original_input_string.len());
        }
        let (_, max_bound) = self.binary_search_max_match(search_string);

        (true, min_bound, max_bound)
    }

    pub fn search_if_match(&self, search_string: &[u8]) -> bool {
        self.binary_search_match(search_string)
    }

}