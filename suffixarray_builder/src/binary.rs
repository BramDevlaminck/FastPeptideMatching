use std::cmp::min;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};

const ONE_GIB: usize = 2usize.pow(30);

/// Trait implemented by structs that are binary serializable
/// In our case this is will be a [i64] since the suffix array is a Vec<i64>
pub trait Serializable {

    /// Serializes self into a vector of bytes
    ///
    /// # Returns
    ///
    /// Returns a vector of bytes
    fn serialize(&self) -> Vec<u8>;
}

impl Serializable for [i64] {
    fn serialize(&self) -> Vec<u8> {
        let mut res = vec![];
        self.iter().for_each(|entry|
            res.extend_from_slice(&entry.to_le_bytes())
        );
        res
    }
}

/// Deserializes a vector of bytes into the suffix array
///
/// # Arguments
/// * `data` - The raw bytes needed to be serialized into a suffix array
///
/// # Returns
///
/// Returns the suffix array, a Vec<i64>
fn deserialize_sa(data: &[u8]) -> Vec<i64> {
    let mut res = vec![];
    if data.len() % 8 != 0 {
        panic!("Serialized data is not a multiple of 8 bytes long!")
    }
    for start in (0..data.len()).step_by(8) {
        res.push(i64::from_le_bytes(data[start..start + 8].try_into().unwrap()));
    }
    res
}

/// Writes the given suffix array with the `sparseness_factor` factor to the given file
///
/// # Arguments
/// * `sparseness_factor` - The sparseness factor of the suffix array
/// * `suffix_array` - The suffix array
/// * `filename` - The name of the file we want to write the suffix array to
///
/// # Returns
///
/// Returns () if writing away the suffix array succeeded
///
/// # Errors
///
/// Returns an io::Error if writing away the suffix array failed
pub fn write_suffix_array(sparseness_factor: u8, suffix_array: &[i64], filename: &str) -> Result<(), std::io::Error> {
    // create the file
    let mut f = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true) // if the file already exists, empty the file
        .open(filename)?;
    f.write_all(&[sparseness_factor])?; // write the sample rate as the first byte

    // write 1 GiB at a time, to minimize extra used memory since we need to translate i64 to [u8; 8]
    let sa_len = suffix_array.len();
    for start_index in (0..sa_len).step_by(ONE_GIB/8) {
        let end_index = min(start_index + ONE_GIB/8, sa_len);
        f.write_all(&suffix_array[start_index..end_index].serialize())?;
    }

    Ok(())
}

/// Loads the suffix array from the file with the given `filename`
///
/// # Arguments
/// * `filename` - The filename of the file where the suffix array is stored
///
/// # Returns
///
/// Returns the sample rate of the suffix array, together with the suffix array
///
/// # Errors
///
/// Returns any error from opening the file or reading the file
pub fn load_suffix_array(filename: &str) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    let mut file = &File::open(filename)?;
    let mut sparseness_factor_buffer = [0_u8; 1];
    file.read_exact(&mut sparseness_factor_buffer).map_err(|_| "Could not read the sample rate from the binary file")?;
    let sparseness_factor = sparseness_factor_buffer[0];

    let mut sa = vec![];
    loop {
        let mut buffer = vec![];
        // use take in combination with read_to_end to ensure that the buffer will be completely filled (except when the file is smaller than the buffer)
        let count = file.take(ONE_GIB as u64).read_to_end(&mut buffer)?;
        if count == 0 {
            break;
        }
        sa.extend_from_slice(&deserialize_sa(&buffer[..count]));
    }

    Ok((sparseness_factor, sa))
}


#[cfg(test)]
mod tests {
    use crate::binary::{deserialize_sa, Serializable};

    #[test]
    fn test_serialize_deserialize() {
        let data: Vec<i64> = vec![5, 2165487362, -12315135];
        let serialized = data.serialize();
        let deserialized = deserialize_sa(serialized.as_ref());
        assert_eq!(data, deserialized);
    }

    #[test]
    fn test_serialize_deserialize_empty() {
        let data: Vec<i64> = vec![];
        let serialized = data.serialize();
        let deserialized = deserialize_sa(serialized.as_ref());
        assert_eq!(data, deserialized);
    }
}