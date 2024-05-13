use std::fs::File;
use std::io;
use std::io::BufRead;
use std::path::Path;
use std::time::{SystemTime, SystemTimeError, UNIX_EPOCH};

use crate::sa_searcher::Searcher;

/// Gets the current time in ms
///
/// # Returns
///
/// Returns the current time in ms
/// 
/// # Errors
/// 
/// Returns a SystemTimeError if getting the current time somehow fails
#[allow(unused)]
pub fn get_time_ms() -> Result<f64, SystemTimeError> {
    Ok(SystemTime::now().duration_since(UNIX_EPOCH)?.as_nanos() as f64 * 1e-6)
}

/// Times how long the function `f`, that has the searcher as only argument executed
///
/// # Returns
///
/// Returns the execution time of `f`in ms
///
/// # Errors
///
/// Returns a SystemTimeError if getting the start or end time failed
#[allow(unused)]
pub fn time_execution(
    searcher: &mut Searcher,
    f: &dyn Fn(&mut Searcher) -> bool,
) -> Result<(bool, f64), SystemTimeError> {
    let start_ms = get_time_ms()?;
    let found = f(searcher);
    let end_ms = get_time_ms()?;
    Ok((found, end_ms - start_ms))
}

/// Opens `filename` and creates an iterator over it per line
///
/// # Arguments
/// * `filename` - The file we want to iterate over per line
/// 
/// # Returns
///
/// Returns an Iterator to the Reader of the lines of the file.
pub fn read_lines<P>(filename: P) -> io::Result<io::Lines<io::BufReader<File>>>
    where P: AsRef<Path>, {
    let file = File::open(filename)?;
    Ok(io::BufReader::new(file).lines())
}
