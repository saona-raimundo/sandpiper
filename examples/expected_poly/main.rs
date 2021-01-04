use average::Variance;
use indicatif::{ProgressBar, ProgressStyle};
use rand::prelude::*;
use rayon::prelude::*;
use sandpiper::{N_REDNECK, N_SANDPIPER, U};
use csv::Writer;
use std::fs::OpenOptions;
use std::fs::File;

// Model parameters
const MUS: [f64; 30] = [-0.0500000000000000, -0.0372751426827093, -0.0277887252403267,  -0.0207165739660756, -0.0154442650096474, -0.0115137436372836,  -0.0085835287378376, -0.0063990451684807, -0.0047705064337644,  -0.0035564261597470, -0.0026513258509018, -0.0019765709878144,  -0.0014735393118657, -0.0010985277619675, -0.0008189555813651,  -0.0006105337229237, -0.0004551546326917, -0.0003393190775256,  -0.0002529633405947, -0.0001885848922832, -0.0001405905753532,  -0.0001048106751227, -0.0000781366573974, -0.0000582511010648,  -0.0000434263620723, -0.0000323744768487, -0.0000241352648763,  -0.0000179929088390, -0.0000134137648850, -0.0000100000000000];
const SIGMAS: [f64; 30] = [0.0000100000000000, 0.0000134137648850, 0.0000179929088390,  0.0000241352648763, 0.0000323744768487, 0.0000434263620723,  0.0000582511010648, 0.0000781366573974, 0.0001048106751227,  0.0001405905753532, 0.0001885848922832, 0.0002529633405947,  0.0003393190775256, 0.0004551546326917, 0.0006105337229237,  0.0008189555813651, 0.0010985277619675, 0.0014735393118657,  0.0019765709878144, 0.0026513258509018, 0.0035564261597470,  0.0047705064337644, 0.0063990451684807, 0.0085835287378376,  0.0115137436372836, 0.0154442650096474, 0.0207165739660756,  0.0277887252403267, 0.0372751426827093, 0.0500000000000000];
const ALPHAS: [f64; 3] = [0., -2., -4.];
const BETAS: [f64; 6] = [0., 10., 50., 100., 1000., 5000.];
const TOTAL: usize = MUS.len() * SIGMAS.len() * ALPHAS.len() * BETAS.len();

// Simulation parameters
const LOWER_S: f64 = -1.;
const UPPER_S: f64 = 1.;
const VARIANCE_SAMPLES: usize = 1_000;
const ERROR_LIMIT: f64 = 1e-6;

fn main() {
    // Computing redneck
    if false {
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
                            let data = [*location, *scale, *shape, *rate, result.mean(), result.error()];
                            save(data, "redneck", counter).unwrap();
                            progress_bar.inc(1);
                        }
                    }
                }
            }
        }
    }

    // Computing sandpiper
    if false {
        let (start, end) = (1, usize::MAX); // Sub-sample
        let mut counter = 0;
        let progress_bar = my_progress_bar(start, end);
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
                            let data = [*location, *scale, *shape, *rate, result.mean(), result.error()];
                            save(data, "sandpiper", counter).unwrap();
                            progress_bar.inc(1);
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
            for counter in 1..=16_200 {
                collect_record(&mut writer, "redneck", counter).unwrap();
            }
            writer.flush().unwrap();
        }
        if true {
            let target_path = "all_sandpiper_poly.csv";
            let output_file = File::create(target_path).unwrap();
            let mut writer = csv::Writer::from_writer(output_file);

            // Getting data
            for counter in 1..=16_200 {
                collect_record(&mut writer, "sandpiper", counter).unwrap();
            }
            writer.flush().unwrap();
        }
    }
}

fn collect_record(writer: &mut csv::Writer<File>, bird: &str, counter: usize) -> Result<(), std::io::Error> {
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
	ProgressBar::new(
        (
        	end.min(TOTAL)
			+ 1 - start
		) as u64
    )
    .with_style(
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