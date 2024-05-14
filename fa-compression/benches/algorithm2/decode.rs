use criterion::black_box;
use fa_compression::algorithm2::{
    decode,
    encode,
    CompressionTable
};

use super::util::generate_annotation;

fn generate_encoded_annotations_and_table(count: usize) -> (Vec<u8>, CompressionTable) {
    let mut random = rand::thread_rng();

    let mut compression_table1 = CompressionTable::new();
    let mut compression_table2 = CompressionTable::new();

    let mut annotations = String::new();
    for _ in 0 .. count {
        let annotation = generate_annotation(&mut random);
        annotations.push_str(&annotation);
        annotations.push(';');
        compression_table1.add_entry(annotation.clone());
        compression_table2.add_entry(annotation);
    }

    annotations.pop();

    (encode(annotations.as_str(), compression_table1), compression_table2)
}

pub fn decode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("decode_algorithm2", |b| {
        b.iter_batched(
            || generate_encoded_annotations_and_table(100),
            |(annotations, ct)| black_box(decode(annotations.as_slice(), ct)),
            criterion::BatchSize::SmallInput
        )
    });
}
