use anyhow::Result;
use preexplorer::prelude::*;
use rand::distributions::Distribution;
use rayon::prelude::*;
use sandpiper::prelude::*;

const SAMPLES: usize = 10_000_000;

fn main() -> Result<()> {
    // Single genetic frequency
    if false {
        let population = N_REDNECK;
        let mutation_rate = U;
        let selection = -0.00001;
        let dominance = 0.5;

        let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
        // let rng = rand::thread_rng();

        let realizations = (0..SAMPLES)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| gen_freq.sample(&mut rand::thread_rng()))
            .collect::<Vec<f64>>();
        // gen_freq.sample_iter(rng).take(SAMPLES);

        // pre::Density::new(&realizations)
        //     .set_title(format!("Genetic frequency. N: {}, U: {}, s: {}, h: {}", population, mutation_rate, selection, dominance))
        //     .plot("testing")
        //     .unwrap();

        let mean: f64 = realizations.iter().sum::<f64>() / SAMPLES as f64;

        println!("Mean: {}", mean);
    }

    // Single Beta random variable
    if false {
        let population = 1;
        let mutation_rate = 1e-4;

        let param = 4. * population as f64 * mutation_rate;
        let beta = sandpiper::Beta::new(param, param).unwrap();
        let rng = rand::thread_rng();

        let realizations = beta.sample_iter(rng).take(SAMPLES).map(|x| {
            // Correcting values
            if (x as f64).is_nan() {
                println!("NAN!!");
                0.
            } else {
                x
            }
        });

        pre::Density::new(realizations)
            .set_title("Beta")
            .plot("testing")
            .unwrap();
    }

    // Various genetic frequencies
    // Estimating statistics
    if true {
        let population = 5_00;
        let mutation_rate = 1.2e-5;
        let selections: Vec<f64> = vec![-1e-0, 1e-0];
        let betas: Vec<f64> = vec![0.];

        for selection in selections {
        	println!("selection: {:?}", selection);
            for beta in &betas {
                let dominance = 1. / (1. + (-beta * selection).exp() as f64);
                let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
                let rng = rand::thread_rng();
                let realizations = gen_freq.sample_iter(rng)
                	.take(SAMPLES);
                let mut histo = quantiles::histogram::Histogram::<f64>::new(
                	vec![0., 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1, 0.5, 1.-1e-1, 1.-1e-2, 1.-1e-3, 1.-1e-4, 1.-1e-5, 1.-1e-6, 1.]
                ).unwrap();
				for realization in realizations {
				    histo.insert(realization);
				}
				for (_bound, count) in histo.iter() {
					println!("{:?}", count);
				}
				// println!("{:#?}", histo);
            }
        }
    }

    // Various genetic frequencies and betas
    // Plotting
    if false {
        let population = 500_000;
        let mutation_rate = 1.2e-8;
        let selections: Vec<f64> = vec![6e-5, 4e-5, 2e-5];
        let betas: Vec<f64> = vec![1e1, 1e2, 1e3, 1e4, 1e5];

        let mut densities_vec = Vec::new();
        for selection in selections {
            for beta in &betas {
                let dominance = 1. / (1. + (-beta * selection).exp() as f64);
                let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
                let rng = rand::thread_rng();
                let realizations = gen_freq.sample_iter(rng).take(SAMPLES);
                densities_vec.push(
                    pre::Density::new(realizations)
                        .set_title(format!("s: {}, beta: {}", selection, beta))
                        .to_owned(),
                )
            }
        }

        pre::Densities::new(densities_vec)
            .set_title("Genetic frequency for red-neck")
            .set_xlabel("frequency")
            .set_ylabel("density value")
            .set_xrange(0., 0.015)
            .plot("genetic_freuqncies")
            .unwrap();
    }

    Ok(())
}
