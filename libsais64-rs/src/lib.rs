// ignore errors because of different style in c code and import the c bindings
#![allow(non_upper_case_globals)]
#![allow(non_camel_case_types)]
#![allow(non_snake_case)]
include!(concat!(env!("OUT_DIR"), "/bindings.rs"));

/// Aid function to return None if the exit code is different from 0, otherwise Some(data)
fn handle_exit_code<T>(exit_code: i64, data: T) -> Option<T> {
    if exit_code == 0 {
        Some(data)
    } else {
        None
    }
}

pub fn sais64(text: &[u8]) -> Option<Vec<i64>> {
    let mut sa = vec![0; text.len()];
    let exit_code = unsafe { libsais64(text.as_ptr(), sa.as_mut_ptr(), text.len() as i64, 0, std::ptr::null_mut()) };
    handle_exit_code(exit_code, sa)
}

pub fn calculate_plcp(text: &[u8], sa: &[i64]) -> Option<Vec<i64>> {
    let mut plcp = vec![0; text.len()];
    let exit_code = unsafe { libsais64_plcp(text.as_ptr(), sa.as_ptr(), plcp.as_mut_ptr(), text.len() as i64) };
    handle_exit_code(exit_code, plcp)
}

pub fn calculate_lcp(text: &[u8], sa: &[i64]) -> Option<Vec<i64>> {
    let plcp = calculate_plcp(text, sa)?;
    let mut lcp = vec![0; text.len()];
    let exit_code = unsafe { libsais64_lcp(plcp.as_ptr(), sa.as_ptr(), lcp.as_mut_ptr(), text.len() as i64) };
    handle_exit_code(exit_code, lcp)
}

#[cfg(test)]
mod tests {
    use crate::{calculate_lcp, calculate_plcp, sais64};

    #[test]
    fn check_build_sa_with_libsais64() {
        let text = "banana$";
        let sa = sais64(text.as_bytes());
        assert_eq!(sa, Some(vec![6, 5, 3, 1, 0, 4, 2]));
    }

    #[test]
    fn check_plcp() {
        let text = "banana$";
        let sa = sais64(text.as_bytes()).expect("Building SA in plcp test failed");
        let plcp = calculate_plcp(text.as_bytes(), &sa);
        assert_eq!(plcp, Some(vec![0, 3, 2, 1, 0, 0, 0]));
    }

    #[test]
    fn check_lcp() {
        let text = "banana$";
        let sa = sais64(text.as_bytes()).expect("Building SA in lcp test failed");
        let lcp = calculate_lcp(text.as_bytes(), &sa);
        assert_eq!(lcp, Some(vec![0, 0, 1, 3, 0, 0, 2]));
    }
}
