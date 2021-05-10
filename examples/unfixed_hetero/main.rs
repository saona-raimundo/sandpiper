use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() -> anyhow::Result<()> {
    // Density of genetic frequency for the UnfixedHeterozygosity
    // Plotting
    if false {
        //Paramteres
        let population = 100;
        let mutation_rate = 1.2e-4;
        let selection: f64 = 5e-1;
        let rate: f64 = 3000.;
        // Random variable
        let hetero = UnfixedHeterozygosity::new(
            population, 
            mutation_rate, 
            Selection::Fixed(selection), 
            Dominance::Sigmoid{rate}, 
            UpperBound::Smallest
        )?;
        // Sampling
        let samples = 1_000;
        let mut rng = rand::thread_rng();
        let realizations: Vec<f64> = (0..samples).map(|_| hetero.sample_frequency(&mut rng)).collect();
        // Reporting
        println!("Empirical max: {:?}", realizations.iter().map(|x| ordered_float::NotNan::new(*x).unwrap()).max());
        println!("Theoretical max: {:?}", 1. - 1. / (2. * population as f64));
        //Plotting
        pre::Density::new(realizations)
            .set_title(format!("s: {}, beta: {}", selection, rate))
            .set_xlabel("frequency")
            .set_ylabel("density value")
            .plot("unfixed_genetic_frequency")
            .unwrap();
    }

    // Fixed parameter computation of expected polymorphisms
    // Printing
    if true {
        // Parameters
        let population_size = 500;
        let mutation_rate = 12e-6;
        let selection = 0.0;
        let beta = 3_000.0;
        // Random variable
        let hetero = UnfixedHeterozygosity::new(
            population_size,
            mutation_rate,
            Selection::Fixed(selection),
            Dominance::Sigmoid{rate: beta}, 
            UpperBound::Smallest
        )?;
        // Computation
        let result = hetero.mc_approx_mean(1000, 1e-4);
        // Reporting
        println!("{:?}", hetero);
        println!("expected Heterozygosity");
        println!("{} \n{}", result.mean(), result.error());
    }

    // Fixed grid parameter computation of expected polymorphisms
    // Data saving
    if false {
        // Parameters
        let population_size = 750;
        let mutation_rate = 1.2e-5;
        let selections = vec![1e-3, 1e-2, 1e-1];
        let betas = vec![0., 1_000.0, 3_000.0];
        // Iterations
        let mut data: Vec<f64> = vec![];
        for selection in selections {
        	println!("{:?}", chrono::offset::Local::now());
        	for beta in &betas {
        		// Random variable
		        let hetero = UnfixedHeterozygosity::new(
		            population_size,
		            mutation_rate,
		            Selection::Fixed(selection),
		            Dominance::Sigmoid{rate: *beta}, 
            		UpperBound::Smallest
		        )?;
		        // Computation
		        let result = hetero.mc_approx_mean(1000, 1e-4);
		        // Collecting
	            data.extend(&[population_size as f64, mutation_rate, *beta, selection, result.mean(), result.error()]);
        	}
        }
        // Save
        // pre::clean().unwrap();
        pre::Data::new(data, 6)
            .save_with_id("expected_polymorphisms_unfixed_allele_frequency")?;
    }

    Ok(())
}