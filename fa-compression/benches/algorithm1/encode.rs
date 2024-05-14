use criterion::black_box;
use fa_compression::algorithm1::encode;

use super::util::generate_annotation;

/// Generate a random number of decoded annotations.
fn generate_decoded_annotations(count: usize) -> String {
    let mut random = rand::thread_rng();

    let mut annotations = String::new();
    for _ in 0 .. count {
        annotations.push_str(&generate_annotation(&mut random));
        annotations.push(';');
    }

    annotations.pop();
    annotations
}

pub fn encode_benchmark(c: &mut criterion::Criterion) {
    c.bench_function("encode_algorithm1", |b| {
        b.iter_batched(
            || generate_decoded_annotations(100),
            |annotations| black_box(encode(annotations.as_str())),
            criterion::BatchSize::SmallInput
        )
    });
}
