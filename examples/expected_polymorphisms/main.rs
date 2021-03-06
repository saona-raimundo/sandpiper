use average::Variance;
use indicatif::{ProgressBar, ProgressStyle};
use preexplorer::prelude::*;
use quadrature::integrate;
use rand::prelude::*;
use rayon::prelude::*;
use sandpiper::{N_REDNECK, N_SANDPIPER, U};

fn main() {
    // Computing the normalizing constant C_{s, \beta}.
    if false {
        // Parameters
        let beta: f64 = 0.;
        let s: f64 = 0.;
        // Computing
        //s >= 0 and N_REDNECK
        let h = |s: f64| 1. / (1. + (-beta * s).exp());
        let integrand = |x: f64| -> f64 {
            (2. * N_REDNECK as f64 * s * (x.powi(2) + 2. * h(s) * x * (1. - x) - 1.)).exp()
                * x.powf(4. * N_REDNECK as f64 * U - 1.)
                * (1. - x).powf(4. * N_REDNECK as f64 * U - 1.)
        };
        let o = integrate(integrand, 0.0, 1.0, 1e-6);
        // Plotting
        let grid = ndarray::Array::linspace(0., 1., 100);
        (&grid, grid.iter().map(|x| integrand(*x)))
            .preexplore()
            .plot("testing")
            .unwrap();
        // Checking

        let reference = 83.25708017910746;
        println!("refernce: {:?}", reference);
        println!(
            "direct integral: {:?}",
            grid.iter().map(|x| integrand(*x)).sum::<f64>() / 100.
        );
        println!("{:?}", o.integral);
        println!("There are problems when evaluation close to the borders!!");
        assert!((o.integral - reference).abs() <= 1e-6);
    }
    // One Monte Carlo approximation
    if false {
        // Parameters
        let location = -4e-5; // mu
        let scale = 1e-5; // sigma
        let shape = -6.; // alpha
        let rate = 1e5; // beta
        let samples = 10000;
        // Computing
        let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
        let result = (0..samples)
            .collect::<Vec<usize>>()
            .into_par_iter()
            .map(|_| {
                let mut rng = thread_rng();
                let s = selection.sample(&mut rng);
                let h = 1. / (1. + (-rate * s).exp());
                let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
                let x = gen_freq.sample(&mut rng);

                // Corrected sampling over machine-presicion
                if x.is_nan() {
                    // Because of the underlying simulation of gamma
                    0.
                } else {
                    2. * x * (1. - x)
                }
            })
            .sum::<f64>()
            / samples as f64;

        println!("Monte Carlo: {:?}", result);
    }
    // Monte Carlo with a given empirical variance
    if false {
        let now = std::time::Instant::now();
        // Parameters
        let location = -0.0000600000; // mu
        let scale = 0.0000100000; // sigma
        let shape = 0.; // alpha
        let rate = 0.; // beta
        let variance_samples = 1000;
        let error_limit = 1e-6; // five significant digits
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
                            let s = selection.sample(&mut rng);
                            let h = 1. / (1. + (-rate * s).exp());
                            let gen_freq =
                                sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
                            let x = gen_freq.sample(&mut rng);

                            // Corrected sampling over machine-presicion
                            if x.is_nan() {
                                // Because of the underlying simulation of gamma
                                0.
                            } else {
                                2. * x * (1. - x)
                            }
                        })
                        .sum::<f64>()
                        / samples as f64;

                    result
                })
                .collect();
            println!("samples: {:?}", samples);
            println!("error: {:?}", variance.error());
        }

        println!("Monte Carlo approximation: {:?}", variance.mean());
        println!("The computations took {} secs.", now.elapsed().as_secs());
    }
    // Approximation for a grid set of parameters
    if false {
        // Parameters
        // mu
        let locations = vec![
            -0.2000000000,
            -0.1523191619,
            -0.1160056354,
            -0.0883494058,
            -0.0672865373,
            -0.0512451448,
            -0.0390280876,
            -0.0297236279,
            -0.0226373905,
            -0.0172405417,
            -0.0131303243,
            -0.0100000000,
        ];
        // sigma
        let scales = vec![1.0];
        // let shapes_and_rates = vec![(0., 1000.), (-2.5, 1000.), (0., 7000.)]; //
        // alpha
        let shapes = vec![0., -2., -4.];
        // beta
        let rates = vec![0., 10., 50., 100., 1000., 5000.]; // vec![1e3]; //
        let variance_samples = 1000;
        let error_limit = 1e-6;
        // Computing redneck
        if true {
            let (start, end) = (1, u64::MAX); //(locations.len() * scales.len() * shapes.len() * rates.len()) as u64); // 2000);
            println!("Redneck");
            let mut data = Vec::new();
            let mut counter = 0;
            let progress_bar = ProgressBar::new(
                u64::min(
                    end,
                    (locations.len() * scales.len() * shapes.len() * rates.len()) as u64,
                ) + 1
                    - start,
            )
            .with_style(
                ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
            );
            for location in locations.iter() {
                for scale in scales.iter() {
                    // for (shape, rate) in shapes_and_rates.iter() {
                    for shape in shapes.iter() {
                        for rate in &rates {
                            counter += 1;
                            if start <= counter && counter <= end {
                                let result: Variance = approximate_expectation_redneck(
                                    *location,
                                    *scale,
                                    *shape,
                                    *rate,
                                    variance_samples,
                                    error_limit,
                                );
                                data.push(*location);
                                data.push(*scale);
                                data.push(*shape);
                                data.push(*rate);
                                data.push(result.mean());
                                data.push(result.error());
                                pre::Data::new(
                                    vec![
                                        *location,
                                        *scale,
                                        *shape,
                                        *rate,
                                        result.mean(),
                                        result.error(),
                                    ],
                                    6,
                                )
                                .set_title(format!(
                                    "Computed value {} of expected polymorphisms. redneck",
                                    counter
                                ))
                                .plot_later(format!("redneck_poly_{}", counter))
                                .unwrap();
                                progress_bar.inc(1);
                            }
                        }
                    }
                }
            }
            pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Redneck")
                .plot_later("redneck_poly")
                .unwrap();
        }

        // Computing sandpiper
        if true {
            let (start, end) = (1, u64::MAX); //(locations.len() * scales.len() * shapes.len() * rates.len()) as u64); // 2000);

            println!("Sandpiper");
            let mut data = Vec::new();
            let mut counter = 0;
            let progress_bar = ProgressBar::new(
                u64::min(
                    end,
                    (locations.len() * scales.len() * shapes.len() * rates.len()) as u64,
                ) + 1
                    - start,
            )
            .with_style(
                ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
            );
            for location in locations.iter() {
                for scale in scales.iter() {
                    // for (shape, rate) in shapes_and_rates.iter() {
                    for shape in shapes.iter() {
                        for rate in &rates {
                            counter += 1;
                            if start <= counter && counter <= end {
                                let result: Variance = approximate_expectation_sandpiper(
                                    *location,
                                    *scale,
                                    *shape,
                                    *rate,
                                    variance_samples,
                                    error_limit,
                                );
                                data.push(*location);
                                data.push(*scale);
                                data.push(*shape);
                                data.push(*rate);
                                data.push(result.mean());
                                data.push(result.error());
                                pre::Data::new(
                                    vec![
                                        *location,
                                        *scale,
                                        *shape,
                                        *rate,
                                        result.mean(),
                                        result.error(),
                                    ],
                                    6,
                                )
                                .set_title(format!(
                                    "Computed value {} of expected polymorphisms. Sandpiper",
                                    counter
                                ))
                                .plot_later(format!("sandpiper_poly_{}", counter))
                                .unwrap();
                                progress_bar.inc(1);
                            }
                        }
                    }
                }
            }
            pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Sandpiper")
                .plot_later("sandpiper_poly")
                .unwrap();
        }
    }
    // Approximation of conditional expectation for a grid set of parameters
    if true {
        // Parameters
        // mu
        let locations = ndarray::Array::geomspace(-0.05, -0.00001, 20).unwrap();
            // vec![];
        // sigma
        let scales = ndarray::Array::geomspace(0.00001, 0.05, 20).unwrap();
            // vec![];
        // alpha
        let shapes = vec![0., -2., -4.];
        // beta
        let rates = vec![0., 100., 1000., 2000., 3000., 5000., 7000.];
        let variance_samples = 1000;
        let error_limit = 1e-4;
        // Computing redneck
        if true {
            let (start, end) = (1, u64::MAX); //(locations.len() * scales.len() * shapes.len() * rates.len()) as u64); // 2000);
            let (lower_bound, upper_bound) = (-1.0, 1.0);
            let mut data = Vec::new();
            let mut counter = 0;
            let progress_bar = ProgressBar::new(
                u64::min(
                    end,
                    (locations.len() * scales.len() * shapes.len() * rates.len()) as u64,
                ) + 1
                    - start,
            )
            .with_style(
                ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
            );
            for location in &locations {
                for scale in &scales {
                    for shape in &shapes {
                        for rate in &rates {
                            counter += 1;
                            if start <= counter && counter <= end {
                                let result: Variance = approximate_conditional_expectation_redneck(
                                    lower_bound,
                                    upper_bound,
                                    *location,
                                    *scale,
                                    *shape,
                                    *rate,
                                    variance_samples,
                                    error_limit,
                                );
                                data.push(*location);
                                data.push(*scale);
                                data.push(*shape);
                                data.push(*rate);
                                data.push(result.mean());
                                data.push(result.error());
                                pre::Data::new(vec![*location, *scale, *shape, *rate, result.mean(), result.error()], 6)
									.set_title(format!("Computed value {} of conditional expected polymorphisms (between {} and {}). redneck", counter, lower_bound, upper_bound))
									.plot_later(format!("redneck_poly_{}", counter))
									.unwrap();
                                progress_bar.inc(1);
                            }
                        }
                    }
                }
            }
            pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Redneck")
                .plot_later("redneck_poly")
                .unwrap();
        }

        // Computing sandpiper
        if true {
            let (start, end) = (1, u64::MAX); //(locations.len() * scales.len() * shapes.len() * rates.len()) as u64); // 2000);
            let (lower_bound, upper_bound) = (-1.0, 1.0);

            let mut data = Vec::new();
            let mut counter = 0;
            let progress_bar = ProgressBar::new(
                u64::min(
                    end,
                    (locations.len() * scales.len() * shapes.len() * rates.len()) as u64,
                ) + 1
                    - start,
            )
            .with_style(
                ProgressStyle::default_bar().template("[{wide_bar}], {pos}/{len} {eta_precise})"),
            );
            for location in &locations {
                for scale in &scales {
                    for shape in &shapes {
                        for rate in &rates {
                            counter += 1;
                            if start <= counter && counter <= end {
                                let result: Variance =
                                    approximate_conditional_expectation_sandpiper(
                                        lower_bound,
                                        upper_bound,
                                        *location,
                                        *scale,
                                        *shape,
                                        *rate,
                                        variance_samples,
                                        error_limit,
                                    );
                                data.push(*location);
                                data.push(*scale);
                                data.push(*shape);
                                data.push(*rate);
                                data.push(result.mean());
                                data.push(result.error());
                                pre::Data::new(vec![*location, *scale, *shape, *rate, result.mean(), result.error()], 6)
									.set_title(format!("Computed value {} of conditional expected polymorphisms (between {} and {}). sandpiper", counter, lower_bound, upper_bound))
									.plot_later(format!("sandpiper_poly_{}", counter))
									.unwrap();
                                progress_bar.inc(1);
                            }
                        }
                    }
                }
            }
            pre::Data::new(data, 6)
                .set_title("Computed values of expected polymorphisms. Sandpiper")
                .plot_later("sandpiper_poly")
                .unwrap();
        }
    }
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

fn approximate_expectation_redneck(
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
                        let s = selection.sample(&mut rng);
                        let h = 1. / (1. + (-rate * s).exp());
                        let gen_freq = sandpiper::GeneticFreq::new(N_REDNECK, U, s, h).unwrap();
                        let x = gen_freq.sample(&mut rng);

                        // Corrected sampling over machine-presicion
                        if x.is_nan() {
                            // Because of the underlying simulation of gamma
                            0.
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

fn approximate_expectation_sandpiper(
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
                        let s = selection.sample(&mut rng);
                        let h = 1. / (1. + (-rate * s).exp());
                        let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
                        let x = gen_freq.sample(&mut rng);

                        // Corrected sampling over machine-presicion
                        if x.is_nan() {
                            // Because of the underlying simulation of gamma
                            0.
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
