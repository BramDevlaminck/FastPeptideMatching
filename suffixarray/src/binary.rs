use checksum_tool::calculate_checksum;
use std::cmp::min;
use std::error::Error;
use std::fs::{File, OpenOptions};
use std::io::{Read, Write};
use crate::index_load_error::IndexLoadError;

const ONE_GIB: usize = 2usize.pow(30);

pub trait Serializable {
    fn serialize(&self) -> Vec<u8>;
}

impl Serializable for [i64] {
    fn serialize(&self) -> Vec<u8> {
        let mut res = vec![];
        self.iter()
            .for_each(|entry| res.extend_from_slice(&entry.to_le_bytes()));
        res
    }
}

fn deserialize_sa(data: &[u8]) -> Vec<i64> {
    let mut res = vec![];
    if data.len() % 8 != 0 {
        panic!("Serialized data is not a multiple of 8 bytes long!")
    }
    for start in (0..data.len()).step_by(8) {
        res.push(i64::from_le_bytes(
            data[start..start + 8].try_into().unwrap(),
        ));
    }
    res
}

struct SuffixArrayBinary {
    database_hash: [u8; 8],
    taxonomy_hash: [u8; 8],
    sample_rate: u8,
    sa: Vec<i64>,
}

// from: https://gist.github.com/taylorsmithgg/ba7b070c0964aa8b86d311ab6f8f5508
// https://dev.to/oliverjumpertz/how-to-write-files-in-rust-m06?comments_sort=top
pub fn write_binary(
    sample_rate: u8,
    suffix_array: &[i64],
    name: &str,
    database_file: &str,
    taxonomy_file: &str,
) -> Result<(), Box<dyn Error>> {
    // create the file
    let mut f = OpenOptions::new()
        .create(true)
        .write(true)
        .truncate(true) // if the file already exists, empty the file
        .open(name.to_owned() + "_sa.bin")?;

    // write checksums
    let database_digest = calculate_checksum(database_file)?;
    f.write_all(database_digest.as_ref())?;
    let taxonomy_digest = calculate_checksum(taxonomy_file)?;
    f.write_all(taxonomy_digest.as_ref())?;

    f.write_all(&[sample_rate])?; // write the sample rate as the first byte

    // write 1 GiB at a time, to minimize extra used memory since we need to translate i64 to [u8; 8]
    let sa_len = suffix_array.len();
    for start_index in (0..sa_len).step_by(ONE_GIB / 8) {
        let end_index = min(start_index + ONE_GIB / 8, sa_len);
        f.write_all(&suffix_array[start_index..end_index].serialize())?;
    }

    Ok(())
}

pub fn load_binary(
    name: &str,
    database_file: &str,
    taxonomy_file: &str,
) -> Result<(u8, Vec<i64>), Box<dyn Error>> {
    let current_database_hash = calculate_checksum(database_file)?;
    let current_taxonomy_hash = calculate_checksum(taxonomy_file)?;
    
    // read the SA and deserialize it into a vec of i64
    let sa_file = File::open(name)?;
    let SuffixArrayBinary { database_hash, taxonomy_hash, sample_rate, sa} = read_sa_file(&sa_file)?;
    
    if database_hash != current_database_hash {
        return Err(Box::new(IndexLoadError::new("The provided database does not match the database used to build the loaded SA")));
    } else if taxonomy_hash != current_taxonomy_hash {
        return Err(Box::new(IndexLoadError::new("The provided taxonomy does not match the taxonomy used to build the loaded SA")));
    }
    
    Ok((sample_rate, sa))
}

/// Aid function to exactly fill the buffer and return an error if it fails
fn fill_buffer(
    mut file: &File,
    buffer: &mut [u8],
    error_message: &str,
) -> Result<(), Box<dyn Error>> {
    file.read_exact(buffer).map_err(|_| error_message)?;
    Ok(())
}

fn read_sa_file(file: &File) -> Result<SuffixArrayBinary, Box<dyn Error>> {
    let mut database_digest_buffer = [0_u8; 8];
    fill_buffer(
        file,
        &mut database_digest_buffer,
        "Could not read the database digest from the binary file",
    )?;
    let mut taxonomy_digest_buffer = [0_u8; 8];
    fill_buffer(
        file,
        &mut taxonomy_digest_buffer,
        "Could not read the taxonomy digest from the binary file",
    )?;

    let mut sample_rate_buffer = [0_u8; 1];
    fill_buffer(
        file,
        &mut sample_rate_buffer,
        "Could not read the sample rate from the binary file",
    )?;
    let sample_rate = sample_rate_buffer[0];

    let mut sa = vec![];
    loop {
        let mut buffer = vec![];
        // use take in combination with read_to_end to ensure that the buffer will be completely filled (except when the file is smaller than the buffer)
        let count = file.take(ONE_GIB as u64).read_to_end(&mut buffer)?;
        if count == 0 {
            break;
        }
        sa.extend(&deserialize_sa(&buffer));
    }

    Ok(SuffixArrayBinary {
        database_hash: database_digest_buffer,
        taxonomy_hash: taxonomy_digest_buffer,
        sample_rate,
        sa,
    })
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
