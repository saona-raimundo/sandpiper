use criterion::{black_box, criterion_group, criterion_main, Criterion};

fn fibonacci(n: u64) -> u64 {
    n + 1
}

pub fn criterion_benchmark(c: &mut Criterion) {
    c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));
}

criterion_group!(benches, criterion_benchmark);
criterion_main!(benches);

// use criterion::{criterion_group, criterion_main, Criterion, BenchmarkId};

// const N: u64 = 500;
// const U: f64 = 1.2e-5;
// const S: f64 = 1e-2;
// const B: f64 = 3000.;

// fn unfixed_hetero(samples: usize) -> Vec<f64> {
//     vec![1.]
//     // let hetero = UnfixedHeterozygosity::new(
//     //     N, 
//     //     U, 
//     //     Selection::Fixed(S), 
//     //     Dominance::Sigmoid{rate: B}
//     // ).unwrap();
//     // // Sampling
//     // let mut rng = rand::thread_rng();
//     // (0..samples).map(|_| hetero.sample_frequency(&mut rng)).collect()
// }

// fn hetero() -> f64 {
//     vec![1.]
//     // let hetero = Heterozygosity::new(
//     //     N, 
//     //     U, 
//     //     Selection::Fixed(S), 
//     //     Dominance::Sigmoid{rate: B}
//     // ).unwrap();
//     // // Sampling
//     // let mut rng = rand::thread_rng();
//     // (0..samples).map(|_| hetero.sample_frequency(&mut rng)).collect()
// }


// fn bench_fibs(c: &mut Criterion) {
//     let mut group = c.benchmark_group("Simulation");
//     for i in [20].iter() {
//         group.bench_with_input(BenchmarkId::new("Unfixed", i), i, 
//             |b, i| b.iter(|| unfixed_hetero(*i)));
//         group.bench_with_input(BenchmarkId::new("Free", i), i, 
//             |b, i| b.iter(|| hetero(*i)));
//     }
//     group.finish();
// }

// criterion_group!(benches, bench_fibs);
// criterion_main!(benches);