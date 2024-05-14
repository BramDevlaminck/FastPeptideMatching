use criterion::criterion_group;

use super::util;

mod decode;
mod encode;

criterion_group!(benches, encode::encode_benchmark, decode::decode_benchmark);
