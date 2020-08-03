use noisy_float::prelude::R64;
use num_traits::float::Float;
use preexplorer::prelude::*;
use rand::distributions::Distribution;
use sandpiper::Parameters;
use sandpiper::Substitutions;
use statrs::statistics::Mean;

const N: u64 = 2_000;

const SAMPLES: u64 = 10_000_000;

fn main() {
    let mu = -1.;
    let sigma = 2.;
    let alpha = 10.;
    let beta = 1000.0;

    let parameters = Parameters::new(mu, sigma, alpha, beta).unwrap();
    let subs = Substitutions::new(N, parameters);

    let mut rng = rand::thread_rng();
    let result: Vec<R64> = (0..SAMPLES)
        .map(|_| subs.sample(&mut rng))
        .filter(|x| x.is_finite())
        .collect();

    pre::Density::new(&result)
        .set_title("Substitutions, empirical distribution.")
        .set_cloud(false)
        .plot("simulation")
        .unwrap();

    let result: Vec<f64> = result.into_iter().map(|x| x.into()).collect();
    println!("Constants");
    println!("N: {}", N);
    println!("Paramenters");
    println!(
        "mu: {}, sigma: {}, alpha: {}, beta: {}",
        mu, sigma, alpha, beta
    );
    println!("Samples: {} \nMean: {}", SAMPLES, result.mean());
}
