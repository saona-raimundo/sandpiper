use crate::constants::*;
use average::Variance;
use csv::Writer;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use rayon::prelude::*;
use sandpiper::{N_REDNECK, N_SANDPIPER, U};
use std::fs::File;
use std::fs::OpenOptions;

pub fn individual_main() {
    // Computing redneck
    if true {
        let (start, end) = (1, usize::MAX); // Sub-sample
        let mut counter: usize = 0;
        let progress_bar = my_progress_bar(start, end);
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
                            let data = [
                                *location,
                                *scale,
                                *shape,
                                *rate,
                                result.mean(),
                                result.error(),
                            ];
                            save(data, "redneck", counter).unwrap();
                            progress_bar.inc(1);
                            progress_bar.reset_eta();
                        }
                    }
                }
            }
        }
    }

    // Computing sandpiper
    if true {
        let (start, end) = (1, usize::MAX); // Sub-sample
        let mut counter = 0;
        let progress_bar = my_progress_bar(start, end);
        for location in &MUS {
            for scale in &SIGMAS {
                for shape in &ALPHAS {
                    for rate in &BETAS {
                        counter += 1;
                        if start <= counter && counter <= end {
                            let result: Variance = approximate_conditional_expectation_sandpiper(
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
                            let data = [
                                *location,
                                *scale,
                                *shape,
                                *rate,
                                result.mean(),
                                result.error(),
                            ];
                            save(data, "sandpiper", counter).unwrap();
                            progress_bar.inc(1);
                            progress_bar.reset_eta();
                        }
                    }
                }
            }
        }
    }

    // Putting all results together
    if true {
        if true {
            let target_path = "all_redneck_poly.csv";
            let output_file = File::create(target_path).unwrap();
            let mut writer = csv::Writer::from_writer(output_file);

            // Getting data
            for counter in 1..=TOTAL {
                collect_record(&mut writer, "redneck", counter).unwrap();
            }
            writer.flush().unwrap();
        }
        if true {
            let target_path = "all_sandpiper_poly.csv";
            let output_file = File::create(target_path).unwrap();
            let mut writer = csv::Writer::from_writer(output_file);

            // Getting data
            for counter in 1..=TOTAL {
                collect_record(&mut writer, "sandpiper", counter).unwrap();
            }
            writer.flush().unwrap();
        }
    }
}

fn collect_record(
    writer: &mut csv::Writer<File>,
    bird: &str,
    counter: usize,
) -> Result<(), std::io::Error> {
    let source_path = format!("{}_poly_{}.csv", bird, counter);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(source_path)?;

    if let Some(result) = rdr.records().next() {
        let record = result?;
        writer.write_record(&record)?;
    }
    Ok(())
}

fn my_progress_bar(start: usize, end: usize) -> ProgressBar {
    ProgressBar::new((end.min(TOTAL) + 1 - start) as u64).with_style(
        ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
    )
}

fn save(data: [f64; 6], bird: &str, counter: usize) -> anyhow::Result<()> {
    let file = OpenOptions::new()
        .append(true)
        .create(true)
        .open(format!("{}_poly_{}.csv", bird, counter))?;
    let mut writer = Writer::from_writer(file);
    writer.serialize(data)?;
    writer.flush()?;
    Ok(())
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
