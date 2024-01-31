use std::cmp::min;
use std::fs::OpenOptions;
use std::io::{Read, Write};

pub trait Serializable {
    fn serialize(&self) -> Vec<u8>;
}

pub trait Deserializable {
    fn deserialize(data: Vec<u8>) -> Vec<i64>;
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

impl Deserializable for Vec<i64> {
    fn deserialize(data: Vec<u8>) -> Vec<i64> {
        let mut res = vec![];
        if data.len() % 8 != 0 {
            panic!("Serialized data is not a multiple of 8 bytes long!")
        }
        for start in (0..data.len()).step_by(8) {
            res.push(i64::from_le_bytes(data[start..start + 8].try_into().unwrap()));
        }
        res
    }
}

// from: https://gist.github.com/taylorsmithgg/ba7b070c0964aa8b86d311ab6f8f5508
// https://dev.to/oliverjumpertz/how-to-write-files-in-rust-m06?comments_sort=top
pub fn write_binary(suffix_array: &Vec<i64>, text: &Vec<u8>, name: &str) -> Result<(), std::io::Error> {
    // create the file
    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(name.to_owned() + "_sa.bin")
        .unwrap();

    // write 1 GiB at a time
    let one_GiB = 1073741824;
    let sa_len = suffix_array.len();
    for start_index in (0..sa_len).step_by(one_GiB) {
        let end_index = min(start_index + one_GiB, sa_len);
        // TODO: remove unwrap and handle error
        f.write_all(&suffix_array[start_index..end_index].serialize()).unwrap();
    }


    let mut f = OpenOptions::new()
        .create(true)
        .append(true)
        .open(name.to_owned() + "_text.bin")
        .unwrap();

    let text_len = text.len();
    for start_index in (0..text_len).step_by(one_GiB) {
        let end_index = min(start_index + one_GiB, text_len);
        // TODO: remove unwrap and handle error
        f.write_all(&text[start_index..end_index]).unwrap();
    }

    Ok(())
}

pub fn load_binary(name: &str) -> Result<(Vec<i64>, Vec<u8>), std::io::Error> {
    let mut f = OpenOptions::new()
        .read(true)
        .open(name.to_owned() + "_sa.bin")
        .unwrap();

    let num_bytes_sa = f.metadata().unwrap().len() as usize;
    let mut buffer_sa = vec![0; num_bytes_sa];

    let read_result = f.read_to_end(&mut buffer_sa).unwrap();
    let sa = Vec::deserialize(buffer_sa);

    let mut f = OpenOptions::new()
        .read(true)
        .open(name.to_owned() + "_text.bin")
        .unwrap();

    let num_bytes_text = f.metadata().unwrap().len() as usize;
    let mut text = vec![0; num_bytes_text];

    let read_result = f.read_to_end(&mut text).unwrap();

    Ok((sa, text))
}


#[cfg(test)]
mod tests {
    use crate::binary::{Deserializable, Serializable};

    #[test]
    fn test_serialize_deserialize() {
        let data: Vec<i64> = vec![5, 2165487362, -12315135];
        let serialized = data.serialize();
        let deserialized = Vec::deserialize(serialized);
        assert_eq!(data, deserialized);
    }

    #[test]
    fn test_serialize_deserialize_empty() {
        let data: Vec<i64> = vec![];
        let serialized = data.serialize();
        let deserialized = Vec::deserialize(serialized);
        assert_eq!(data, deserialized);
    }
}