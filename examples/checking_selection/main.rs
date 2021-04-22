use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() {
    // Fixed parameter computation
    if false {
        // Parameters
        let population_size = 5_000;
        let mutation_rate = 1.2e-6;
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
    if false {
        // Parameters
        let population_size = 5_000;
        let mutation_rate = 1.2e-6;
        let dominance = 0.5; 
        let selections = [1e-1]; // ndarray::Array::geomspace(1e-7, 1e-2, 6).unwrap();
        // Iterate
        let mut data: Vec<f64> = Vec::new();
        for selection in selections.iter() {
            // Compute
            let hetero = Heterozygosity::new(
                    population_size,
                    mutation_rate,
                    Selection::Fixed(*selection),
                    Dominance::Fixed(dominance),
                ).expect("Could not construct Heterozygosity");
            let result: average::Variance = hetero.mc_approx_mean(1000, 1e-5);
            // Recover
            data.extend(&[population_size as f64, mutation_rate, dominance, *selection, result.mean(), result.error()]);
        }
        // Save
        pre::Data::new(data, 6)
            .set_title("Computed values of expected polymorphisms. Redneck")
            .save_with_id("varying_selection")
            .unwrap();
    }
}