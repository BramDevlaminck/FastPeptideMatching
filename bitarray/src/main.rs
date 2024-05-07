use std::io::{stdin, stdout, BufWriter, Read, Write};

use bitarray::{binary::Binary, BitArray};

pub fn main() {
    let stdin = stdin();
    let stdout = stdout();

    let mut reader = stdin.lock().bytes();
    //let mut writer = BufWriter::new(stdout.lock());

    let sample_rate = reader.next().unwrap().unwrap();
    //writer.write(&[sample_rate]).unwrap();

    // eprintln!("Reading the sample rate from the binary file: {}", sample_rate);

    //let size: usize = 29_401_012_218;
    let size: usize = 25_000;
    let mut bitarray = BitArray::<37>::with_capacity(size);
    //writer.write(&size.to_le_bytes()).unwrap();

    let mut index = 0;
    let mut buffer = [0; 8];
    let mut byte_count = 0;
    for byte in reader {
        buffer[byte_count] = byte.unwrap();
        byte_count += 1;

        if byte_count == 8 {
            //println!("{:064b}", i64::from_le_bytes(buffer));
            bitarray.set(index, u64::from_le_bytes(buffer));
            byte_count = 0;
            index += 1;
        }
    }

    // eprintln!("value (1023): {}", bitarray.get(1023));
    // eprintln!("value (1024): {}", bitarray.get(1024));
    // eprintln!("value (1025): {}", bitarray.get(1025));
    // eprintln!("value (1999): {}", bitarray.get(1999));
    // //eprintln!("value ({}): {}", amount_of_entries - 1, sa.get(amount_of_entries as usize - 1));

    // for i in (1_000 .. 5_000).step_by(250) {
    //     eprintln!("value ({}): {}", i, bitarray.get(i as usize));
    // }

    //bitarray.write_binary(&mut writer).unwrap();
}
