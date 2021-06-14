use average::Variance;
use constants::*;
use csv::Writer;
use indicatif::{ProgressBar, ProgressStyle};
use sandpiper::prelude::*;
use std::fs::File;
use std::fs::OpenOptions;

mod constants {
    // Model parameters
    pub const MUS: [f64; 30] = [
        -0.0500000000000000,
        -0.0372751426827093,
        -0.0277887252403267,
        -0.0207165739660756,
        -0.0154442650096474,
        -0.0115137436372836,
        -0.0085835287378376,
        -0.0063990451684807,
        -0.0047705064337644,
        -0.0035564261597470,
        -0.0026513258509018,
        -0.0019765709878144,
        -0.0014735393118657,
        -0.0010985277619675,
        -0.0008189555813651,
        -0.0006105337229237,
        -0.0004551546326917,
        -0.0003393190775256,
        -0.0002529633405947,
        -0.0001885848922832,
        -0.0001405905753532,
        -0.0001048106751227,
        -0.0000781366573974,
        -0.0000582511010648,
        -0.0000434263620723,
        -0.0000323744768487,
        -0.0000241352648763,
        -0.0000179929088390,
        -0.0000134137648850,
        -0.0000100000000000,
    ];
    pub const SIGMAS: [f64; 30] = [
        0.0000100000000000,
        0.0000134137648850,
        0.0000179929088390,
        0.0000241352648763,
        0.0000323744768487,
        0.0000434263620723,
        0.0000582511010648,
        0.0000781366573974,
        0.0001048106751227,
        0.0001405905753532,
        0.0001885848922832,
        0.0002529633405947,
        0.0003393190775256,
        0.0004551546326917,
        0.0006105337229237,
        0.0008189555813651,
        0.0010985277619675,
        0.0014735393118657,
        0.0019765709878144,
        0.0026513258509018,
        0.0035564261597470,
        0.0047705064337644,
        0.0063990451684807,
        0.0085835287378376,
        0.0115137436372836,
        0.0154442650096474,
        0.0207165739660756,
        0.0277887252403267,
        0.0372751426827093,
        0.0500000000000000,
    ];
    pub const ALPHAS: [f64; 3] = [0., -2., -4.];
    pub const BETAS: [f64; 5] = [0., 1000., 3000., 5000., 7000.];
    pub const TOTAL: usize = MUS.len() * SIGMAS.len() * ALPHAS.len() * BETAS.len();

    // Simulation parameters
    pub const UPPER_GEN_FREQ: sandpiper::UpperBound = sandpiper::UpperBound::Smallest;
    pub const VARIANCE_SAMPLES: usize = 1_000;
    pub const ERROR_LIMIT: f64 = 1e-6;
}

fn main() -> anyhow::Result<()> {
    println!("Started!");
    let redneck_bool = true;
    let sandpiper_bool = true;

    simulate(redneck_bool, sandpiper_bool)?;
    gather_records(redneck_bool, sandpiper_bool)?;
    Ok(())
}

fn simulate(redneck_bool: bool, sandpiper_bool: bool) -> anyhow::Result<()> {
    let args: Vec<String> = std::env::args().collect();
    let start: usize = args[1].parse().unwrap();
    let end = start + 1000 - 1;
    println!("Computing from {} to {}", start, end);
    let mut counter: usize = 0;
    let progress_bar = my_progress_bar(start, end);
    for location in &MUS {
        for scale in &SIGMAS {
            for shape in &ALPHAS {
                for rate in &BETAS {
                    counter += 1;
                    if start <= counter && counter <= end {
                        if redneck_bool {
                            let result: Variance = approximate_conditional_expectation(
                                sandpiper::N_REDNECK,
                                sandpiper::U,
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
                            save(data, "redneck", counter)?;
                        }
                        if sandpiper_bool {
                            let result: Variance = approximate_conditional_expectation(
                                sandpiper::N_SANDPIPER,
                                sandpiper::U,
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
                            save(data, "sandpiper", counter)?;
                        }

                        // Report progress
                        println!(
                            "{} of {}. Done in {} hours. ETA: {} hours",
                            counter,
                            end.min(TOTAL) - start + 1,
                            progress_bar.elapsed().as_secs() / 3600,
                            progress_bar.eta().as_secs() / 3600
                        );
                        println!("");
                        progress_bar.inc(1);
                    }
                }
            }
        }
    }

    Ok(())
}

fn gather_records(redneck_bool: bool, sandpiper_bool: bool) -> anyhow::Result<()> {
    if redneck_bool {
        let target_path = "all_redneck_poly.csv";
        let output_file = File::create(target_path)?;
        let mut writer = csv::Writer::from_writer(output_file);

        // Getting data
        for counter in 1..=TOTAL {
            collect_record(&mut writer, "redneck", counter)?;
        }
        writer.flush().unwrap();
    }
    if sandpiper_bool {
        let target_path = "all_sandpiper_poly.csv";
        let output_file = File::create(target_path).unwrap();
        let mut writer = csv::Writer::from_writer(output_file);

        // Getting data
        for counter in 1..=TOTAL {
            collect_record(&mut writer, "sandpiper", counter)?;
        }
        writer.flush()?;
    }
    Ok(())
}

//////////////////////////////////////////////////////////////

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

fn approximate_conditional_expectation(
    population_size: u64,
    mutation_rate: f64,
    location: f64,
    scale: f64,
    shape: f64,
    rate: f64,
    variance_samples: usize,
    error_limit: f64,
) -> Variance {
    let hetero = UnfixedHeterozygosity::new(
        population_size,
        mutation_rate,
        Selection::SkewNormal {
            location,
            scale,
            shape,
            bounds: Some((-1., 1.)),
        },
        Dominance::Sigmoid { rate },
        UPPER_GEN_FREQ,
    )
    .unwrap();

    hetero.mc_approx_mean(variance_samples, error_limit)
}
