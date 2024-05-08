use std::io::{stdin, stdout, BufWriter, Read, Write};

use bitarray::{binary::Binary, BitArray};

pub fn main() {
    let stdin = stdin();
    let stdout = stdout();

    let mut reader = stdin.lock();
    let mut writer = BufWriter::new(stdout.lock());

    let mut sample_rate_buffer = [0_u8; 1];
    reader.read_exact(&mut sample_rate_buffer).map_err(|_| "Could not read the sample rate from the binary file").unwrap();
    writer.write(&sample_rate_buffer).unwrap();

    // eprintln!("Reading the sample rate from the binary file: {}", sample_rate);

    //let size: usize = 29_401_012_218;
    let size: usize = 200_000_000;
    let mut bitarray = BitArray::<37>::with_capacity(size);
    writer.write(&size.to_le_bytes()).unwrap();

    let mut index = 0;
    let mut buffer = vec![0; 8 * 1024];
    loop {
        let (finished, bytes_read) = fill_buffer(&mut reader, &mut buffer);
        for buffer_slice in buffer[..bytes_read].chunks_exact(8) {
            bitarray.set(index, u64::from_le_bytes(buffer_slice.try_into().unwrap()));
            index += 1;
        }

        if finished {
            break;
        }
    }

    for i in (0 .. 20_000).step_by(250) {
        eprintln!("value ({}): {}", i, bitarray.get(i as usize));
    }

    bitarray.write_binary(&mut writer).unwrap();
}

fn fill_buffer<T: Read>(input: &mut T, buffer: &mut Vec<u8>) -> (bool, usize) {
    // Store the buffer size in advance, because rust will complain
    // about the buffer being borrowed mutably while it's borrowed
    let buffer_size = buffer.len();

    let mut writable_buffer_space = buffer.as_mut();

    loop {
        match input.read(writable_buffer_space) {
            // No bytes written, which means we've completely filled the buffer
            // or we've reached the end of the file
            Ok(0) => {
                return (
                    !writable_buffer_space.is_empty(),
                    buffer_size - writable_buffer_space.len()
                );
            }

            // We've read {bytes_read} bytes
            Ok(bytes_read) => {
                // Shrink the writable buffer slice
                writable_buffer_space = writable_buffer_space[bytes_read..].as_mut();
            }

            Err(err) => {
                panic!("Error while reading input: {}", err);
            }
        }
    }
}
