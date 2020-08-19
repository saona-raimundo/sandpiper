use indicatif::{ProgressStyle, ProgressBar};
use quadrature::integrate;
use sandpiper::{N_REDNECK, N_SANDPIPER, U};
use preexplorer::prelude::*;
use rand::prelude::*;
use rayon::prelude::*;
use average::Variance;

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
		    (2. * N_REDNECK as f64 * s * ( x.powi(2) + 2. * h(s) * x * (1. - x) - 1.)).exp() 
		    * x.powf(4. * N_REDNECK as f64 * U - 1.)
		    * (1. - x).powf(4. * N_REDNECK as f64 * U - 1.) 
		};
		let o = integrate(integrand , 0.0, 1.0, 1e-6);
		// Plotting
		let grid = ndarray::Array::linspace(0., 1., 100);
		(&grid, grid.iter().map(|x| integrand(*x))).preexplore()
			.plot("testing")
			.unwrap();
		// Checking

		let reference = 83.25708017910746;
		println!("refernce: {:?}", reference);
		println!("direct integral: {:?}", grid.iter().map(|x| integrand(*x)).sum::<f64>() / 100.);
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
		let result = (0..samples).collect::<Vec<usize>>()
			.into_par_iter()
			.map(|_| {
				let mut rng = thread_rng();
				let s = selection.sample(&mut rng);
				let h = 1. / (1. + (-rate * s).exp());
				let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
				let x = gen_freq.sample(&mut rng);

				// Corrected sampling over machine-presicion
				if x.is_nan() { // Because of the underlying simulation of gamma
					0. 
				} else {
					2. * x * (1. - x)
				}
			})
			.sum::<f64>() / samples as f64;

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
			variance = (0..variance_samples).map(|_| {
				// Computing
				let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
				let result = (0..samples).collect::<Vec<usize>>()
					.into_par_iter()
					.map(|_| {
						let mut rng = thread_rng();
						let s = selection.sample(&mut rng);
						let h = 1. / (1. + (-rate * s).exp());
						let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
						let x = gen_freq.sample(&mut rng);

						// Corrected sampling over machine-presicion
						if x.is_nan() { // Because of the underlying simulation of gamma
							0. 
						} else {
							2. * x * (1. - x)
						}
					})
					.sum::<f64>() / samples as f64;

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
	if true {	
		// Parameters
		let locations = vec![-0.0080000000, -0.0058363249, -0.0042578361, -0.0031062643, -0.0022661460, -0.0016532455, -0.0012061098, -0.0008799061, -0.0006419272, -0.0004683120, -0.0003416526, -0.0002492494, -0.0001818376, -0.0001326579, -0.0000967793, -0.0000706045, -0.0000515088, -0.0000375778, -0.0000274145, -0.0000200000]; // vec![-0.0000600000, -0.0000566291, -0.0000534476, -0.0000504448, -0.0000476107, -0.0000449359, -0.0000424113, -0.0000400285, -0.0000377797, -0.0000356571, -0.0000336539, -0.0000317631, -0.0000299786, -0.0000282944, -0.0000267047, -0.0000252044, -0.0000237884, -0.0000224519, -0.0000211905, -0.0000200000]; // vec![-4e-5]; //
		let scales = vec![0.0000100000, 0.0000142165, 0.0000202110, 0.0000287331, 0.0000408486, 0.0000580726, 0.0000825592, 0.0001173706, 0.0001668606, 0.0002372181, 0.0003372423, 0.0004794423, 0.0006816015, 0.0009690021, 0.0013775867, 0.0019584530, 0.0027842444, 0.0039582349, 0.0056272444, 0.0080000000]; // vec![0.0000100000, 0.0000117078, 0.0000137073, 0.0000160482, 0.0000187889, 0.0000219977, 0.0000257544, 0.0000301527, 0.0000353022, 0.0000413311, 0.0000483897, 0.0000566536, 0.0000663290, 0.0000776566, 0.0000909188, 0.0001064459, 0.0001246247, 0.0001459081, 0.0001708263, 0.0002000000]; // vec![1e-5]; // 
		let shapes_and_rates = vec![(0., 1000.), (-2.5, 1000.), (0., 7000.)]; // let shapes = vec![0., -2.5, -5., -7.5, -10.]; // vec![-5.]; // 
		// let rates = vec![7000.]; // vec![0., 1000., 3000., 5000.]; // vec![1e3]; // 
		let variance_samples = 1000;
		let error_limit = 1e-6;
		// Computing redneck
		if true {
			let mut data = Vec::new();
			let mut counter = 0;
			let (start, end) = (1, (locations.len()*scales.len()*shapes_and_rates.len()) as u64); // 2000);
			let progress_bar = ProgressBar::new(end + 1 - start).with_style(
				ProgressStyle::default_bar()
					.template("[{wide_bar}], {pos}/{len} {eta_precise})"));
			for location in locations.iter() {
				for scale in scales.iter() {
					for (shape, rate) in shapes_and_rates.iter() {
					// for shape in shapes.iter() {
					// 	for rate in &rates {
							counter += 1;
							if start <= counter && counter <= end {
								let result: Variance = approximate_expectation_redneck(*location, *scale, *shape, *rate, variance_samples, error_limit);
								data.push(*location);
								data.push(*scale);
								data.push(*shape);
								data.push(*rate);
								data.push(result.mean());
								data.push(result.error());
								pre::Data::new(vec![*location, *scale, *shape, *rate, result.mean(), result.error()], 6)
									.set_title(format!("Computed value {} of expected polymorphisms. redneck", counter))
									.plot_later(format!("redneck_poly_{}", counter))
									.unwrap();
								progress_bar.inc(1);

							// }
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
			let mut data = Vec::new();
			let mut counter = 0;
			let (start, end) = (1, (locations.len()*scales.len()*shapes_and_rates.len()) as u64); // 2000);
			let progress_bar = ProgressBar::new(end + 1 - start).with_style(
				ProgressStyle::default_bar()
					.template("[{wide_bar}], {pos}/{len} {eta_precise})"));
			for location in locations.iter() {
				for scale in scales.iter() {
					for (shape, rate) in shapes_and_rates.iter() {
					// for shape in shapes.iter() {
					// 	for rate in &rates {
							counter += 1;
							if start <= counter && counter <= end {
								let result: Variance = approximate_expectation_sandpiper(*location, *scale, *shape, *rate, variance_samples, error_limit);
								data.push(*location);
								data.push(*scale);
								data.push(*shape);
								data.push(*rate);
								data.push(result.mean());
								data.push(result.error());
								pre::Data::new(vec![*location, *scale, *shape, *rate, result.mean(), result.error()], 6)
									.set_title(format!("Computed value {} of expected polymorphisms. Sandpiper", counter))
									.plot_later(format!("sandpiper_poly_{}", counter))
									.unwrap();
								progress_bar.inc(1);
							// }
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
	if false {	
		// Parameters
		let locations = vec![-0.0000600000, -0.0000566291, -0.0000534476, -0.0000504448, -0.0000476107, -0.0000449359, -0.0000424113, -0.0000400285, -0.0000377797, -0.0000356571, -0.0000336539, -0.0000317631, -0.0000299786, -0.0000282944, -0.0000267047, -0.0000252044, -0.0000237884, -0.0000224519, -0.0000211905, -0.0000200000]; // vec![-4e-5]; //
		let scales = vec![0.0000100000, 0.0000117078, 0.0000137073, 0.0000160482, 0.0000187889, 0.0000219977, 0.0000257544, 0.0000301527, 0.0000353022, 0.0000413311, 0.0000483897, 0.0000566536, 0.0000663290, 0.0000776566, 0.0000909188, 0.0001064459, 0.0001246247, 0.0001459081, 0.0001708263, 0.0002000000]; // vec![1e-5]; // 
		let shapes = vec![-10.]; // vec![0., -2.5, -5., -7.5]; // vec![-5.]; // 
		let rates = vec![0., 1000., 3000., 5000.]; // vec![1e3]; // 
		let variance_samples = 1000;
		let error_limit = 1e-6;
		// Computing redneck
		if false {
			let lower_bound = -0.01;
			let upper_bound = 0.0005;
			let mut data = Vec::new();
			let mut counter = 0;
			let (start, end) = (1550, 1600);
			for location in &locations {
				for scale in &scales {
					for shape in &shapes {
						for rate in &rates {
							let now = std::time::Instant::now();	
							counter += 1;
							if start <= counter && counter <= end {
								let result: Variance = approximate_conditional_expectation_redneck(lower_bound, upper_bound, *location, *scale, *shape, *rate, variance_samples, error_limit);
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
								println!("{:?}", counter);
								println!("Remaining {} hours.", (end - counter) * now.elapsed().as_secs() / 3600);	
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
		if false {
			let lower_bound = -0.01;
			let upper_bound = 0.0005;
			let mut data = Vec::new();
			let mut counter = 0;
			let (start, end) = (1520, 1600);
			for location in &locations {
				for scale in &scales {
					for shape in &shapes {
						for rate in &rates {
							let now = std::time::Instant::now();	
							counter += 1;
							if start <= counter && counter <= end {
								let result: Variance = approximate_conditional_expectation_sandpiper(lower_bound, upper_bound, *location, *scale, *shape, *rate, variance_samples, error_limit);
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
								println!("{:?}", counter);
								println!("Remaining {} hours.", (end - counter) * now.elapsed().as_secs() / 3600);	
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

fn approximate_conditional_expectation_sandpiper(lower_bound:f64, upper_bound: f64, location: f64, scale: f64, shape: f64, rate: f64, variance_samples: usize, error_limit: f64) -> Variance {
	// Variance recursion
	let mut variance: Variance = [1., 0., -1.].iter().collect();
	let mut samples = 1000;
	while variance.error() > error_limit {
		samples *= 2;
		variance = (0..variance_samples).map(|_| {
			// Computing
			let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
			let result = (0..samples).collect::<Vec<usize>>()
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
					if x.is_nan() { // Because of the underlying simulation of gamma
						panic!("NaN!");
					} else {
						2. * x * (1. - x)
					}
				})
				.sum::<f64>() / samples as f64;

			result
			})
			.collect();
	}
	variance
}


fn approximate_conditional_expectation_redneck(lower_bound:f64, upper_bound: f64, location: f64, scale: f64, shape: f64, rate: f64, variance_samples: usize, error_limit: f64) -> Variance {
	// Variance recursion
	let mut variance: Variance = [1., 0., -1.].iter().collect();
	let mut samples = 1000;
	while variance.error() > error_limit {
		samples *= 2;
		variance = (0..variance_samples).map(|_| {
			// Computing
			let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
			let result = (0..samples).collect::<Vec<usize>>()
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
					if x.is_nan() { // Because of the underlying simulation of gamma
						panic!("NaN!");
					} else {
						2. * x * (1. - x)
					}
				})
				.sum::<f64>() / samples as f64;

			result
			})
			.collect();
	}
	variance
}

fn approximate_expectation_redneck(location: f64, scale: f64, shape: f64, rate: f64, variance_samples: usize, error_limit: f64) -> Variance {
	// Variance recursion
	let mut variance: Variance = [1., 0., -1.].iter().collect();
	let mut samples = 1000;
	while variance.error() > error_limit {
		samples *= 2;
		variance = (0..variance_samples).map(|_| {
			// Computing
			let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
			let result = (0..samples).collect::<Vec<usize>>()
				.into_par_iter()
				.map(|_| {
					let mut rng = thread_rng();
					let s = selection.sample(&mut rng);
					let h = 1. / (1. + (-rate * s).exp());
					let gen_freq = sandpiper::GeneticFreq::new(N_REDNECK, U, s, h).unwrap();
					let x = gen_freq.sample(&mut rng);

					// Corrected sampling over machine-presicion
					if x.is_nan() { // Because of the underlying simulation of gamma
						0. 
					} else {
						2. * x * (1. - x)
					}
				})
				.sum::<f64>() / samples as f64;

			result
			})
			.collect();
	}
	variance
}

fn approximate_expectation_sandpiper(location: f64, scale: f64, shape: f64, rate: f64, variance_samples: usize, error_limit: f64) -> Variance {
	// Variance recursion
	let mut variance: Variance = [1., 0., -1.].iter().collect();
	let mut samples = 1000;
	while variance.error() > error_limit {
		samples *= 2;
		variance = (0..variance_samples).map(|_| {
			// Computing
			let selection = sandpiper::SkewNormal::new(location, scale, shape).unwrap();
			let result = (0..samples).collect::<Vec<usize>>()
				.into_par_iter()
				.map(|_| {
					let mut rng = thread_rng();
					let s = selection.sample(&mut rng);
					let h = 1. / (1. + (-rate * s).exp());
					let gen_freq = sandpiper::GeneticFreq::new(N_SANDPIPER, U, s, h).unwrap();
					let x = gen_freq.sample(&mut rng);

					// Corrected sampling over machine-presicion
					if x.is_nan() { // Because of the underlying simulation of gamma
						0. 
					} else {
						2. * x * (1. - x)
					}
				})
				.sum::<f64>() / samples as f64;

			result
			})
			.collect();
	}
	variance
}