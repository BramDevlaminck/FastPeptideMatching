use std::io::{stdin, stdout, BufReader, BufWriter, Read};

use bitarray::{binary::Binary, BitArray};

pub fn main() {
    let stdin = stdin();
    let stdout = stdout();

    let mut reader = BufReader::new(stdin.lock());
    let mut writer = BufWriter::new(stdout.lock());

    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file").unwrap();

    let mut bitarray = BitArray::<38>::with_capacity(30_000_000_000);

    // this buffer is 1GiB big
    let mut index = 0;
    let mut buffer = [0; 8 * 1024];
    let mut bytes_read = reader.read(&mut buffer).unwrap();
    while bytes_read > 0 {
        for buffer_slice in buffer.chunks_exact(8) {
            bitarray.set(index, u64::from_le_bytes(buffer_slice.try_into().unwrap()));
            index += 1;
        }
        bytes_read = reader.read(&mut buffer).unwrap();
    }

    bitarray.write_binary(&mut writer).unwrap();
}
