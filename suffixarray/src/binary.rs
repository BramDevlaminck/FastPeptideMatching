use std::fs::OpenOptions;
use std::io::{Read, Write};


pub trait Serializable {

    fn serialize(&self) -> Vec<u8>;

    fn deserialize(data: Vec<u8>) -> Self;

}

impl Serializable for Vec<i64> {
    fn serialize(&self) -> Vec<u8> {
        let mut res = vec![];
        self.iter().for_each(|entry|
            res.extend_from_slice(&entry.to_le_bytes())
        );
        res
    }

    fn deserialize(data: Vec<u8>) -> Self {
        let mut res = vec![];
        if data.len() % 8 != 0 {
            panic!("Serialized data is not a multiple of 8 bytes long!")
        }
        for start in (0..data.len()).step_by(8) {
            res.push(i64::from_le_bytes(data[start..start+8].try_into().unwrap()));
        }
        res
    }

}

// from: https://gist.github.com/taylorsmithgg/ba7b070c0964aa8b86d311ab6f8f5508
pub fn write_binary(suffix_array: &Vec<i64>, text: &[u8], name: &str) -> Result<(), std::io::Error> {
    let mut f = OpenOptions::new()
        .create(true)
        .write(true)
        .read(true)
        .open(name.to_owned() + "_sa.bin")
        .unwrap();

    let write_result = f.write(&suffix_array.serialize()).unwrap();

    let mut f = OpenOptions::new()
        .create(true)
        .write(true)
        .read(true)
        .open(name.to_owned() + "_text.bin")
        .unwrap();

    let write_result = f.write(text).unwrap();

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
    use crate::binary::Serializable;

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