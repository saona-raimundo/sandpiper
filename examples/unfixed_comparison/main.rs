use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() -> anyhow::Result<()> {
    // Fixed parameter computation of expected polymorphisms
    // Printing
    if false {
        // Parameters
        let population_size = 750;
        let mutation_rate = 8e-6;
        let selection = 0.00666666;
        let beta = 0.0; // 3_000.0;
        let upper_bound = UpperBound::Smallest;
        // Random variable
        let hetero = UnfixedHeterozygosity::new(
            population_size,
            mutation_rate,
            Selection::Fixed(selection),
            Dominance::Sigmoid{rate: beta}, 
            upper_bound
        )?;
        // Computation
        let result = hetero.mc_approx_mean(1000, 1e-4);
        // Reporting
        println!("{:?}", hetero);
        println!("expected Heterozygosity");
        println!("{} \n{}", result.mean(), result.error());
    }

    // Fixed parameter comparison of expected polymorphisms (varying upper_bound)
    // Printing
    if true {
        // Parameters
        let population_size = 500;
        let mutation_rate = 12e-6;
        let selection = 0.0;
        let beta = 0.; // 3_000.0;
        // Reporting
        println!("N = {}, U = {}, s = {}, beta = {}", population_size, mutation_rate, selection, beta);
        println!("expected Heterozygosity");
        // Iteration
        let upper_bounds = vec![UpperBound::Smallest, UpperBound::Midpoint, UpperBound::Largest];
        for upper_bound in upper_bounds {
            // Random variable
            let hetero = UnfixedHeterozygosity::new(
                population_size,
                mutation_rate,
                Selection::Fixed(selection),
                Dominance::Sigmoid{rate: beta}, 
                upper_bound
            )?;
            // Computation
            let result = hetero.mc_approx_mean(1000, 1e-4);
            // Reporting
            println!("{:?}", upper_bound);
            println!("{} \n{}", result.mean(), result.error());        
        }
    
    }

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