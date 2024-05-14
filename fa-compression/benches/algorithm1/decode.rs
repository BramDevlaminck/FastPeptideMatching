use criterion::black_box;
use fa_compression::algorithm1::{
    decode,
    encode
};

use super::util::generate_annotation;

/// Generate a random number of encoded annotations.
fn generate_encoded_annotations(count: usize) -> Vec<u8> {
    let mut random = rand::thread_rng();

    let mut annotations = String::new();
    for _ in 0 .. count {
        annotations.push_str(&generate_annotation(&mut random));
        annotations.push(';');
    }
    annotations.pop();

    encode(annotations.as_str())
}

pub fn decode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("decode_algorithm1", |b| {
        b.iter_batched(
            || generate_encoded_annotations(100),
            |annotations| black_box(decode(annotations.as_slice())),
            criterion::BatchSize::SmallInput
        )
    });
}
