use anyhow::Result;
use preexplorer::prelude::*;
use rand::distributions::Distribution;
use rayon::prelude::*;
use sandpiper::prelude::*;

const SAMPLES: usize = 10_000; // _000;

fn main() -> Result<()> {
    // Fixed parameter genetic frequency
    // Printing mean
    if false {
        // Parameters
        let population = N_REDNECK;
        let mutation_rate = U;
        let selection = -0.00001;
        let dominance = 0.5;
        // Random variable
        let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
        // Simulation
        let realizations = (0..SAMPLES)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| gen_freq.sample(&mut rand::thread_rng()))
            .collect::<Vec<f64>>();
        // Reporting
        let mean: f64 = realizations.iter().sum::<f64>() / SAMPLES as f64;
        println!("Mean: {}", mean);
    }

    // Fixed parameter genetic frequency
    // Plotting and saving histogram
    if true {
        // Parameters
        let population = 5_000;
        let mutation_rate = 1.2e-6;
        let selection = 0.000_1;
        let beta: f64 = 0.; //3_000.;
                            // Random variable
        let dominance = 1. / (1. + (-beta * selection).exp());
        let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
        // Simulation
        let samples = 100_000_000;
        let realizations = (0..samples)
            .collect::<Vec<usize>>()
            .par_iter()
            .map(|_| gen_freq.sample(&mut rand::thread_rng()))
            .collect::<Vec<f64>>();
        // Reporting
        if false {
            // Plotting
            pre::Density::new(realizations.clone())
                .set_title(format!(
                    "Genertic frequency, N = {}, U = {}, s = {}, beta = {}",
                    population, mutation_rate, selection, beta
                ))
                .set_cdf(false)
                .set_pdf(false)
                .plot(format!(
                    "genertic frequency, N = {}, U = {}, s = {}, beta = {}",
                    population, mutation_rate, selection, beta
                ))
                .unwrap();
        }
        if true {
            // Histogram
            let grid = vec![
                0.,
                1e-6,
                1e-5,
                0.5e-4,
                1e-4,
                1e-3,
                1e-2,
                1e-1,
                0.5,
                1. - 1e-1,
                1. - 1e-2,
                1. - 1e-3,
                1. - 1e-4,
                1. - 0.5e-4,
                1. - 1e-5,
                1. - 1e-6,
                1.,
            ];
            // Computations
            let init_samples = 1000;
            let repetitions = 1000;
            let error = 1e-4;
            let histo = sandpiper::distribution::helper::approx_histogram(
                gen_freq,
                &mut rand::thread_rng(),
                grid.clone(),
                init_samples,
                repetitions,
                error,
            )
            .unwrap();
            // Reporting
            let mut data = vec![];
            for i in 0..grid.len() {
                data.push(grid[i]);
                data.push(
                    histo.total_below(quantiles::histogram::Bound::Finite(grid[i])) as f64
                        / histo.count() as f64,
                );
            }
            let title = format!(
                "density frequency, beta = {}, N = {}, U = {}, s = {}",
                beta, population, mutation_rate, selection
            );
            println!("{:?}", title);
            println!("{:#?}", data);
            pre::Data::new(data, 2)
                .set_title(&title)
                .save_with_id(&title)
                .unwrap();
        }
    }

    // Fixed Beta random variable
    // Plotting density
    if false {
        // Parameters
        let population = 1;
        let mutation_rate = 1e-4;
        // Random variable
        let param = 4. * population as f64 * mutation_rate;
        let beta = sandpiper::Beta::new(param, param).unwrap();
        let rng = rand::thread_rng();
        // Simulation
        let realizations = beta.sample_iter(rng).take(SAMPLES);
        //Reporting
        pre::Density::new(realizations)
            .set_title("Beta")
            .plot("testing")
            .unwrap();
    }

    // Various genetic frequencies
    // Printing estimating statistics
    if false {
        // Parameters
        let population = 5_00;
        let mutation_rate = 1.2e-5;
        let selections: Vec<f64> = vec![-1e-0, 1e-0];
        let betas: Vec<f64> = vec![0.];
        // Simualtion
        for selection in selections {
            println!("selection: {:?}", selection);
            for beta in &betas {
                let dominance = 1. / (1. + (-beta * selection).exp() as f64);
                let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance)?;
                let rng = rand::thread_rng();
                // Simulation
                let realizations = gen_freq.sample_iter(rng).take(SAMPLES);
                // Reporting
                let mut histo = quantiles::histogram::Histogram::<f64>::new(vec![
                    0.,
                    1e-6,
                    1e-5,
                    1e-4,
                    1e-3,
                    1e-2,
                    1e-1,
                    0.5,
                    1. - 1e-1,
                    1. - 1e-2,
                    1. - 1e-3,
                    1. - 1e-4,
                    1. - 1e-5,
                    1. - 1e-6,
                    1.,
                ])
                .unwrap();
                for realization in realizations {
                    histo.insert(realization);
                }
                for (_bound, count) in histo.iter() {
                    println!("{:?}", count);
                }
            }
        }
    }

    // Various genetic frequencies and betas
    // Plotting
    if false {
        // Plotting
        let population = 500;
        let mutation_rate = 1.2e-4;
        let selections: Vec<f64> = vec![1e-6, 1e-3, 1e-1, 5e-1];
        let betas: Vec<f64> = vec![0., 1e1, 1e5];
        // Simulation
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
        // Reporting
        pre::Densities::new(densities_vec)
            .set_title("Genetic frequency for red-neck")
            .set_xlabel("frequency")
            .set_ylabel("density value")
            .plot("genetic_freuqncies")
            .unwrap();
    }

    // One genetic frequency for the UnfixedHeterozygosity
    // Plotting
    if false {
        // Parameters
        let population = 100;
        let mutation_rate = 1.2e-4;
        let selection: f64 = 5e-1;
        let beta: f64 = 0.;
        // Random variable
        let dominance = 1. / (1. + (-beta * selection).exp() as f64);
        let hetero = UnfixedHeterozygosity::new(
            population,
            mutation_rate,
            sandpiper::Selection::Fixed(selection),
            sandpiper::Dominance::Fixed(dominance),
            sandpiper::UpperBound::Smallest,
        )?;
        // Simulation
        let realizations = (0..SAMPLES).map(|_| hetero.sample_frequency(&mut rand::thread_rng()));
        // Reporting
        pre::Density::new(realizations)
            .set_title(format!("s: {}, beta: {}", selection, beta))
            .set_xlabel("frequency")
            .set_ylabel("density value")
            .plot("unfixed_genetic_frequency")
            .unwrap();
    }

    // Various genetic frequencies and betas for the UnfixedHeterozygosity
    // Plotting
    if false {
        // Parameters
        let population = 500;
        let mutation_rate = 1.2e-4;
        let selections: Vec<f64> = vec![1e-6, 1e-3, 1e-1, 5e-1];
        let betas: Vec<f64> = vec![0., 1e1, 1e5];
        // Simulation
        let mut densities_vec = Vec::new();
        for selection in selections {
            for beta in &betas {
                let dominance = 1. / (1. + (-beta * selection).exp() as f64);
                let hetero = UnfixedHeterozygosity::new(
                    population,
                    mutation_rate,
                    sandpiper::Selection::Fixed(selection),
                    sandpiper::Dominance::Fixed(dominance),
                    sandpiper::UpperBound::Smallest,
                )?;
                let mut rng = rand::thread_rng();
                let realizations = (0..SAMPLES).map(|_| hetero.sample_frequency(&mut rng));
                densities_vec.push(
                    pre::Density::new(realizations)
                        .set_title(format!("s: {}, beta: {}", selection, beta))
                        .to_owned(),
                )
            }
        }
        // Reporting
        pre::Densities::new(densities_vec)
            .set_title("Genetic frequency for red-neck")
            .set_xlabel("frequency")
            .set_ylabel("density value")
            .plot("unfixed_genetic_frequencies")
            .unwrap();
    }

    Ok(())
}
