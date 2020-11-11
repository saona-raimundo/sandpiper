use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() {
    // Fixed parameter computation
    if false {
        let hetero = Heterozygosity::new(
            N_REDNECK,
            U,
            Selection::Fixed(1e-7),
            // SkewNormal {
            //     location: 0.00001,
            //     scale: 0.0001,
            //     shape: -2.,
            //     bounds: None,
            // },
            Dominance::Fixed(0.5),
        )
        .unwrap();

        println!("{:?}", hetero.mc_approx_mean(1000, 1e-6));
    }

    // Varying s for fixed parameter
    if true {
        let selections = 
        	ndarray::Array::linspace(-1., 1., 20);
        	// [-1e2, -5e1, -1e1, -5e0, -2e0, -15e-1, -10e-1, -5e-1, -1e-1, -1e-2, 1e-2, 1e-1, 5e-1, 10e-1, 15e-1, 2e0, 5e0, 1e1, 5e1, 1e2];
        let values: Vec<f64> = selections.iter().map(|selection| {
            let hetero = Heterozygosity::new(
                N_REDNECK,
                U,
                Selection::Fixed(*selection),
                Dominance::Sigmoid{rate: 50.},
            )
            .unwrap();

            hetero.mc_approx_mean(1000, 1e-6).mean()
        })
        .collect();

        (&selections, values)
            .preexplore()
            .set_title("Expected polymorphisms.")
            .set_xlabel("selection, s")
            .plot("varying s")
            .unwrap();
    }
}
