use criterion::black_box;
use fa_compression::algorithm2::{
    encode,
    CompressionTable
};

use super::util::generate_annotation;

fn generate_decoded_annotations_and_table(count: usize) -> (String, CompressionTable) {
    let mut random = rand::thread_rng();

    let mut compression_table = CompressionTable::new();

    let mut annotations = String::new();
    for _ in 0 .. count {
        let annotation = generate_annotation(&mut random);
        annotations.push_str(&annotation);
        annotations.push(';');
        compression_table.add_entry(annotation);
    }

    annotations.pop();

    (annotations, compression_table)
}

pub fn encode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("encode_algorithm2", |b| {
        b.iter_batched(
            || generate_decoded_annotations_and_table(100),
            |(annotations, ct)| black_box(encode(annotations.as_str(), ct)),
            criterion::BatchSize::SmallInput
        )
    });
}
