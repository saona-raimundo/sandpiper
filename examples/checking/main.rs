use preexplorer::prelude::*;
use sandpiper::prelude::*;

fn main() {
    // Fixed parameter computation
    if false {
        // Parameters
        let population_size = 5000;
        let mutation_rate = 1.2e-6;
        let beta = 3_000.0;
        // Random variable
        let hetero = Heterozygosity::new(
            population_size,
            mutation_rate,
            Selection::SkewNormal {
                location: -0.13458588,
                scale: 0.0377358,
                shape: 0.0,
                bounds: Some((-1., 1.)),
            },
            Dominance::Sigmoid { rate: beta },
        )
        .unwrap();
        // Computation
        let result = hetero.mc_approx_mean(1000, 1e-6);
        println!("{:?}", result);
        println!("{} \n{}", result.mean(), result.error());
    }

    // Fixed Dominance for various s
    if false {
        let selections = grid(); // ndarray::Array::linspace(-1e-1, 1e-1, 20);
        let h = Dominance::Sigmoid { rate: 2500.0 };

        // Compute
        let values: Vec<f64> = selections
            .iter()
            .map(|&selection| {
                let hetero = Heterozygosity::new(N_REDNECK, U, Selection::Fixed(selection), h)
                    .expect("Could not construct Heterozygosity");
                hetero.mc_approx_mean(1000, 1e-4).mean()
            })
            .collect();

        // Visualize
        // (&selections, values)
        //     .preexplore()
        //     .set_title(format!("Expected polymorphisms. {}", h))
        //     .set_style(3)
        //     .set_logy(10)
        //     .set_logx(10)
        //     .set_xlabel("selection, s")
        //     .plot(format!("varying s with {}", h))
        //     .unwrap();

        let process = (&selections, values)
            .preexplore()
            .set_title(format!("h = {}", h))
            .set_style(3)
            .to_owned();
        let reference = (
            &selections,
            selections
                .iter()
                .map(|_| EMPIRICAL_MEAN_POLYMORPHISMS_REDNECK),
        )
            .preexplore()
            .set_title("redneck value")
            .to_owned();

        (reference + process)
            .set_title("Expected polymorphisms")
            .set_logy(10)
            .set_logx(10)
            .set_xlabel("selection, s")
            .plot(format!("varying s with {}", h))
            .unwrap();
    }

    // Varying fixed Dominance for various s
    if false {
        let hs = vec![0., 1.];
        let mut all = Vec::new();
        for h in hs {
            // Compute
            let values: Vec<f64> = grid()
                .iter()
                .map(|&selection| {
                    let hetero = Heterozygosity::new(
                        N_REDNECK,
                        U,
                        Selection::Fixed(selection),
                        sandpiper::Dominance::Fixed(h),
                    )
                    .expect("Could not construct Heterozygosity");

                    hetero.mc_approx_mean(1000, 1e-5).mean()
                })
                .collect();
            // Visualize
            all.push(
                (grid().to_vec(), values)
                    .preexplore()
                    .set_title(format!("h = {:?}", h))
                    .to_owned(),
            )
            // .plot(format!("{}", h))
            // .unwrap();
        }
        pre::Processes::new(all)
            .set_title("Expected polymorphisms")
            .set_style(3)
            .set_logy(10)
            .set_logx(10)
            .set_xlabel("selection, s")
            .plot("fixed dominance")
            .unwrap();
    }

    // Sigmoidal plots
    if false {
        let betas = vec![1, 100, 5_000, 100_000];
        for beta in betas {
            // Compute
            let values: Vec<f64> = grid()
                .iter()
                .map(|&selection| 1. / (1. + (-(beta as f64) * selection).exp()))
                .collect();
            // Visualize
            (grid().to_vec(), values)
                .preexplore()
                .set_title(format!("Beta = {}", beta))
                .set_rangey(-0.05, 1.05)
                .set_style(3)
                .set_logx(10)
                .plot(format!("Beta {}", beta))
                .unwrap();
        }
    }

    // Varying sigmoid Dominance for various s
    if false {
        let betas = vec![1, 100, 5_000, 100_000];
        for beta in betas {
            // Compute
            let values: Vec<f64> = grid()
                .iter()
                .map(|&selection| {
                    let hetero = Heterozygosity::new(
                        N_REDNECK,
                        U,
                        Selection::Fixed(selection),
                        sandpiper::Dominance::Sigmoid {
                            rate: (beta as f64),
                        },
                    )
                    .expect("Could not construct Heterozygosity");

                    hetero.mc_approx_mean(1000, 1e-5).mean()
                })
                .collect();
            // Visualize
            (grid().to_vec(), values)
                .preexplore()
                .set_title(format!("Beta = {}", beta))
                .set_rangey(1e-8, 1e-1)
                .set_style(3)
                .set_logx(10)
                .set_logy(10)
                .plot(format!("Sigmoidal {}", beta))
                .unwrap();
        }
    }
}

fn grid() -> ndarray::Array1<f64> {
    ndarray::Array::geomspace(1e-7, 1e-0, 8).unwrap()
}
