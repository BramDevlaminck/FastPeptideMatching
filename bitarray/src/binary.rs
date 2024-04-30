use std::io::{BufRead, Result, Write};

use crate::BitArray;

pub trait Binary {
    fn write_binary<W: Write>(&self, writer: W) -> Result<()>;
    fn read_binary<R: BufRead>(&mut self, reader: R) -> Result<()>;
}

impl<const B: usize> Binary for BitArray<B> {
    fn write_binary<W: Write>(&self, mut writer: W) -> Result<()> {
        for value in self.data.iter() {
            writer.write_all(&value.to_le_bytes())?;
        }

        Ok(())
    }

    fn read_binary<R: BufRead>(&mut self, mut reader: R) -> Result<()> {
        self.data.clear();

        let mut i = 0;

        let mut buffer = [0; 8 * 1024];
        let mut bytes_read = reader.read(&mut buffer)?;
        while bytes_read > 0 {
            for buffer_slice in buffer.chunks_exact(8) {
                let number = u64::from_le_bytes(buffer_slice.try_into().unwrap());
                eprintln!("Read number: {}", number);
                self.data.push(u64::from_le_bytes(buffer_slice.try_into().unwrap()));

                if i == 5 {
                    eprintln!("0: {}, 1: {}, 2: {}", self.get(0), self.get(1), self.get(2));
                    break;
                }
                i+=1;
            }
            break;
            bytes_read = reader.read(&mut buffer)?;
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_write_binary() {
        let mut bitarray = BitArray::<40>::with_capacity(4);
        bitarray.set(0, 0x1234567890);
        bitarray.set(1, 0xabcdef0123);
        bitarray.set(2, 0x4567890abc);
        bitarray.set(3, 0xdef0123456);

        let mut buffer = Vec::new();
        bitarray.write_binary(&mut buffer).unwrap();

        assert_eq!(buffer, vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12,
            0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01,
            0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
        ]);
    }

    #[test]
    fn test_read_binary() {
        let buffer = vec![
            0xef, 0xcd, 0xab, 0x90, 0x78, 0x56, 0x34, 0x12,
            0xde, 0xbc, 0x0a, 0x89, 0x67, 0x45, 0x23, 0x01,
            0x00, 0x00, 0x00, 0x00, 0x56, 0x34, 0x12, 0xf0
        ];

        let mut bitarray = BitArray::<40>::with_capacity(4);
        bitarray.read_binary(&buffer[..]).unwrap();

        assert_eq!(bitarray.get(0), 0x1234567890);
        assert_eq!(bitarray.get(1), 0xabcdef0123);
        assert_eq!(bitarray.get(2), 0x4567890abc);
        assert_eq!(bitarray.get(3), 0xdef0123456);
    }
}
