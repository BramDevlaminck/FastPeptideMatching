use std::error::Error;
use std::fs::File;
use std::io::{Read};

const ONE_GIB: usize = 2usize.pow(30);

pub fn calculate_checksum(file_name: &str) -> Result<[u8; 8], Box<dyn Error>> {
    let file = File::open(file_name)?;
    let mut hash = 0;

    loop {
        let mut buffer = vec![];
        let count = (&file).take(ONE_GIB as u64).read_to_end(&mut buffer)?;
        if count == 0 {
            break;
        }
        // for the last iteration we have to append some zero's to the end of the buffer to ensure that the buffer is always a multiple of 8 bytes long
        let to_append = buffer.len() % 8;
        if to_append != 0 {
            buffer.extend(vec![0; 8 - to_append]);
        }
        // update hash value with current buffer
        buffer.chunks(8).map(|chunk| u64::from_le_bytes(chunk.try_into().unwrap())).for_each(|val| hash ^= val);
    }
    Ok(hash.to_le_bytes())
}
