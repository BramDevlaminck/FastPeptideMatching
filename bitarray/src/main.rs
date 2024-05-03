use std::io::{stdin, stdout, BufReader, BufWriter, Read, Write};

use bitarray::{binary::Binary, BitArray};

pub fn main() {
    let stdin = stdin();
    let stdout = stdout();

    eprintln!("size of usize: {}", std::mem::size_of::<usize>());

    let mut reader = BufReader::new(stdin.lock());
    let mut writer = BufWriter::new(stdout.lock());

    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file").unwrap();
    writer.write(&sample_rate_buffer).unwrap();

    let size: usize = 30_000_000_000;
    let mut bitarray = BitArray::<37>::with_capacity(size);
    writer.write(&size.to_le_bytes()).unwrap();

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
