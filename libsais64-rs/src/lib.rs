// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));


/// Builds the suffix array over the `text` using the libsais64 algorithm
///
/// # Arguments
/// * `text` - The text used for suffix array construction
///
/// # Returns
///
/// Returns Some with the suffix array build over the text if construction succeeds
/// Returns None if construction of the suffix array failed
pub fn sais64(text: &[u8]) -> Option<Vec<i64>> {
    let mut sa = vec![0; text.len()];
    let exit_code = unsafe { libsais64(text.as_ptr(), sa.as_mut_ptr(), text.len() as i64, 0, std::ptr::null_mut()) };
    if exit_code == 0 {
        Some(sa)
    } else {
        None
    }
}

#[cfg(test)]
mod tests {
    use crate::{sais64};

    #[test]
    fn check_build_sa_with_libsais64() {
        let text = "banana$";
        let sa = sais64(text.as_bytes());
        assert_eq!(sa, Some(vec![6, 5, 3, 1, 0, 4, 2]));
    }
}
