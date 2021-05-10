use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() {
    // Fixed parameter computation
    // Fixed dominance
    if false {
        // Parameters
        let population_size = 500;
        let mutation_rate = 1.2e-5;
        let dominance = 0.5;
        let selection = 1e-7;
        // Random variable
        let hetero = Heterozygosity::new(
            population_size,
            mutation_rate,
            Selection::Fixed(selection),
            Dominance::Fixed(dominance),
        )
        .unwrap();
        // Computation
        println!("{:?}", hetero.mc_approx_mean(1000, 1e-6));
    }

    // Fixed parameter computation
    // Sigmoidal dominance
    if false {
        // Parameters
        let population_size = 5000;
        let mutation_rate = 1.2e-6;
        let beta = 3_000.0;
        let selection = 1e-7;
        // Random variable
        let hetero = Heterozygosity::new(
            population_size,
            mutation_rate,
            Selection::Fixed(selection),
            Dominance::Sigmoid{rate: beta},
        )
        .unwrap();
        // Computation
        println!("{:?}", hetero.mc_approx_mean(1000, 1e-6));
    }

    // Varying only selection s
    // Ploting mean polymorphisms
    if false {
        // Parameters
        let population_size = 5_000;
        let mutation_rate = 1.2e-6;
        let dominance = 0.5; 
        let selections = ndarray::Array::geomspace(1e-7, 1e-2, 6).unwrap(); // grid(); // ndarray::Array::linspace(-1e-1, 1e-1, 20);

        // Compute
        let values: Vec<f64> = selections
            .iter()
            .map(|&selection| {
            	println!("{:?}", chrono::offset::Local::now());
                let hetero = Heterozygosity::new(
                    population_size,
                    mutation_rate,
                    Selection::Fixed(selection),
                    Dominance::Fixed(dominance),
                ).expect("Could not construct Heterozygosity");
                hetero.mc_approx_mean(1000, 1e-5).mean()
            })
            .collect();

        // Visualize
        (&selections, values)
            .preexplore()
            .set_title(format!("Expected polymorphisms. h = {}", dominance))
            .set_style(3)
            .set_logy(10)
            .set_logx(10)
            .set_xlabel("selection, s")
            .plot(format!("varying s with {}", dominance))
            .unwrap();
    }

    // Varying only selection s
    // Saving data
    if true {
        // Parameters
        let population_size = 500;
        let mutation_rate = 1.2e-5;
        let rate = 3_000.0; 
        let selections = [0.]; // [-1e-1, -1e-2, -1e-3, -1e-4, -1e-5, -1e-6, -1e-7, 0., 1e-7, 1e-6, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1]; // ndarray::Array::geomspace(1e-7, 1e-2, 6).unwrap();
        // Iterate
        let mut data: Vec<f64> = Vec::new();
        for selection in selections.iter() {
        	println!("{:?}", chrono::offset::Local::now());
            // Compute
            let hetero = Heterozygosity::new(
                    population_size,
                    mutation_rate,
                    Selection::Fixed(*selection),
                    Dominance::Sigmoid{rate},
                ).expect("Could not construct Heterozygosity");
            let result: average::Variance = hetero.mc_approx_mean(1_000, 1e-6);
            // Recover
            data.extend(&[population_size as f64, mutation_rate, rate, *selection, result.mean(), result.error()]);
        }
        // Report
        println!("{:#?}", data);
        // Save
        pre::clean().unwrap();
        pre::Data::new(data, 6)
            .set_title("Computed values of expected polymorphisms. Redneck")
            .save_with_id("varying_selection")
            .unwrap();
    }
}