use criterion::criterion_main;

mod algorithm1;
mod algorithm2;
mod util;

criterion_main!(algorithm1::benches, algorithm2::benches);
