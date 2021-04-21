use preexplorer::prelude::*;
use average::Variance;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use rayon::prelude::*;
use sandpiper::{N_REDNECK, N_SANDPIPER, U};
// const U: f64 = 0.0000012;
// const N_REDNECK: u64 = 5000;
// const N_SANDPIPER: u64 = 5000;

// use csv::Writer;
// use std::fs::File;

use crate::constants::*;

pub fn collective_main() {
    // Computing redneck
    if false {
        let (start, end) = (1, usize::MAX); // Sub-sample
        let mut counter: usize = 0;
        let progress_bar = my_progress_bar(end.min(TOTAL) + 1 - start);
        let mut data: Vec<f64> = Vec::with_capacity(end.min(TOTAL) + 1 - start);
        for location in &MUS {
            for scale in &SIGMAS {
                for shape in &ALPHAS {
                    for rate in &BETAS {
                        counter += 1;
                        if start <= counter && counter <= end {
                            let result: Variance = approximate_conditional_expectation_redneck(
                                LOWER_S,
                                UPPER_S,
                                *location,
                                *scale,
                                *shape,
                                *rate,
                                VARIANCE_SAMPLES,
                                ERROR_LIMIT,
                            );
                            // Save
                            data.extend(&[*location, *scale, *shape, *rate, result.mean(), result.error()]);
                            progress_bar.inc(1);
                        }
                    }
                }
            }
        }
        pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Redneck")
                .plot_later("redneck_all")
                .unwrap();
    }

    // Computing sandpiper
    if true {
        let (start, end) = (1, usize::MAX); // Sub-sample
        let mut counter = 0;
        let progress_bar = my_progress_bar(end.min(TOTAL) + 1 - start);
        let mut data: Vec<f64> = Vec::with_capacity(end.min(TOTAL) + 1 - start);
        for location in &MUS {
            for scale in &SIGMAS {
                for shape in &ALPHAS {
                    for rate in &BETAS {
                        counter += 1;
                        if start <= counter && counter <= end {
                            let result: Variance =
                                approximate_conditional_expectation_sandpiper(
                                    LOWER_S,
                                    UPPER_S,
                                    *location,
                                    *scale,
                                    *shape,
                                    *rate,
                                    VARIANCE_SAMPLES,
                                    ERROR_LIMIT,
                                );
                            // Save
                            data.extend(&[*location, *scale, *shape, *rate, result.mean(), result.error()]);
                            progress_bar.inc(1);
                        }
                    }
                }
            }
        }
        pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Sandpiper")
                .plot_later("sandpiper_all")
                .unwrap();
    }
}

fn my_progress_bar(end: usize) -> ProgressBar {
	ProgressBar::new(end as u64).with_style(
        ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
    )
}

fn approximate_conditional_expectation_sandpiper(
    lower_bound: f64,
    upper_bound: f64,
    location: f64,
    scale: f64,
    shape: f64,
    rate: f64,
    variance_samples: usize,
    error_limit: f64,
) -> Variance {
    // Variance recursion
    let mut variance: Variance = [1., 0., -1.].iter().collect();
    let mut samples = 1000;
    while variance.error() > error_limit {
        samples *= 2;
        variance = (0..variance_samples)
            .map(|_| {
                // Computing
                let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
                let result = (0..samples)
                    .collect::<Vec<usize>>()
                    .into_par_iter()
                    .map(|_| {
                        let mut rng = thread_rng();
                        let mut s = std::f64::INFINITY;
                        while s > upper_bound || s < lower_bound {
                            s = selection.sample(&mut rng);
                        }
                        let h = 1. / (1. + (-rate * s).exp());
                        let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
                        let x = gen_freq.sample(&mut rng);

                        // Corrected sampling over machine-presicion
                        if x.is_nan() {
                            // Because of the underlying simulation of gamma
                            panic!("NaN!");
                        } else {
                            2. * x * (1. - x)
                        }
                    })
                    .sum::<f64>()
                    / samples as f64;

                result
            })
            .collect();
    }
    variance
}

fn approximate_conditional_expectation_redneck(
    lower_bound: f64,
    upper_bound: f64,
    location: f64,
    scale: f64,
    shape: f64,
    rate: f64,
    variance_samples: usize,
    error_limit: f64,
) -> Variance {
    // Variance recursion
    let mut variance: Variance = [1., 0., -1.].iter().collect();
    let mut samples = 1000;
    while variance.error() > error_limit {
        samples *= 2;
        variance = (0..variance_samples)
            .map(|_| {
                // Computing
                let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
                let result = (0..samples)
                    .collect::<Vec<usize>>()
                    .into_par_iter()
                    .map(|_| {
                        let mut rng = thread_rng();
                        let mut s = std::f64::INFINITY;
                        while s > upper_bound || s < lower_bound {
                            s = selection.sample(&mut rng);
                        }
                        let h = 1. / (1. + (-rate * s).exp());
                        let gen_freq = sandpiper::GeneticFreq::new(N_REDNECK, U, s, h).unwrap();
                        let x = gen_freq.sample(&mut rng);

                        // Corrected sampling over machine-presicion
                        if x.is_nan() {
                            // Because of the underlying simulation of gamma
                            panic!("NaN!");
                        } else {
                            2. * x * (1. - x)
                        }
                    })
                    .sum::<f64>()
                    / samples as f64;

                result
            })
            .collect();
    }
    variance
}
