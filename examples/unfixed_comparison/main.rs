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
            Dominance::Sigmoid { rate: beta },
            upper_bound,
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
        let population_size = 5000;
        let mutation_rate = 1.2e-6;
        let selection = 1e-5;
        let beta = 0.; // 3_000.0;
                       // Reporting
        println!(
            "N = {}, U = {}, s = {}, beta = {}",
            population_size, mutation_rate, selection, beta
        );
        println!("expected Heterozygosity");
        // Iteration
        let upper_bounds = vec![
            UpperBound::Smallest,
            UpperBound::Midpoint,
            UpperBound::Largest,
        ];
        for upper_bound in upper_bounds {
            // Random variable
            let hetero = UnfixedHeterozygosity::new(
                population_size,
                mutation_rate,
                Selection::Fixed(selection),
                Dominance::Sigmoid { rate: beta },
                upper_bound,
            )?;
            // Computation
            let result = hetero.mc_approx_mean(1000, 1e-5);
            // Reporting
            println!("{:?}", upper_bound);
            println!("{} \n{}", result.mean(), result.error());
        }
    }

    Ok(())
}
