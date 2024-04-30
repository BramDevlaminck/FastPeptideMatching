pub mod binary;

pub struct BitArray<const B: usize> {
    pub data: Vec<u64>,
    pub mask: u64,
    pub len: usize,
}

impl<const B: usize> BitArray<B> {
    pub fn with_capacity(capacity: usize) -> Self {
        Self {
            data: vec![0; capacity * B / 64 + 1],
            mask: (1 << B) - 1,
            len: capacity,
        }
    }

    pub fn get(&self, index: usize) -> u64 {
        let start_block = index * B / 64;
        let start_block_offset = index * B % 64;

        if index >= 1023 && index <= 1025 {
            eprintln!("start_block: {}, start_block_offset: {}", start_block, start_block_offset);
        }

        if start_block_offset + B <= 64 {
            if index >= 1023 && index <= 1025 {
                eprintln!("data: {0:64b}", self.data[start_block]);
                eprintln!("shift: {}", (64 - start_block_offset - B));
                eprintln!("shifted: {0:64b}", self.data[start_block] >> (64 - start_block_offset - B));
                eprintln!("result: {0:64b}", self.data[start_block] >> (64 - start_block_offset - B) & self.mask);
            }

            return self.data[start_block] >> (64 - start_block_offset - B) & self.mask;
        }

        let end_block = (index + 1) * B / 64;
        let end_block_offset = (index + 1) * B % 64;

        if index >= 1023 && index <= 1025 {
            eprintln!("end_block: {}, end_block_offset: {}", end_block, end_block_offset);
        }

        let a = self.data[start_block] << end_block_offset;
        let b = self.data[end_block] >> (64 - end_block_offset);

        if index >= 1023 && index <= 1025 {
            eprintln!("a: {}, b: {}", a, b);
            eprintln!("result: {}", (a | b) & self.mask);
        }

        (a | b) & self.mask
    }

    pub fn set(&mut self, index: usize, value: u64) {
        let start_block = index * B / 64;
        let start_block_offset = index * B % 64;

        if start_block_offset + B <= 64 {
            self.data[start_block] &= !(self.mask << (64 - start_block_offset - B));
            self.data[start_block] |= value << (64 - start_block_offset - B);
            return;
        }

        let end_block = (index + 1) * B / 64;
        let end_block_offset = (index + 1) * B % 64;

        self.data[start_block] &= !(self.mask >> start_block_offset);
        self.data[start_block] |= value >> end_block_offset;

        self.data[end_block] &= !(self.mask << (64 - end_block_offset));
        self.data[end_block] |= value << (64 - end_block_offset);
    }

    pub fn len(&self) -> usize {
        self.len
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_bitarray_get() {
        let mut bitarray = BitArray::<40>::with_capacity(4);
        bitarray.data = vec![ 0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144eb9ca32eb4a4 ];

        assert_eq!(bitarray.get(0), 0b0001110011111010110001000111111100110010);
        assert_eq!(bitarray.get(1), 0b1100001001010010011000010100110111001001);
        assert_eq!(bitarray.get(2), 0b1111001101001101101101101011101001010001);
        assert_eq!(bitarray.get(3), 0b0000100010010001010001001110101110011100);
    }

    #[test]
    fn test_bitarray_set() {
        let mut bitarray = BitArray::<40>::with_capacity(4);
        bitarray.data = vec![ 0, 0, 0 ];

        bitarray.set(0, 0b0001110011111010110001000111111100110010);
        bitarray.set(1, 0b1100001001010010011000010100110111001001);
        bitarray.set(2, 0b1111001101001101101101101011101001010001);
        bitarray.set(3, 0b0000100010010001010001001110101110011100);

        assert_eq!(bitarray.data, vec![ 0x1cfac47f32c25261, 0x4dc9f34db6ba5108, 0x9144EB9C00000000 ]);
    }
}
