//! Visualization of fitting in parameter space for mean substitutions.

use itertools::Itertools;
use num_traits::real::Real;
use preexplorer::prelude::*;
use rayon::prelude::*;
use sandpiper::{prelude::*, ExpBinary};
use statrs::statistics::Mean;

fn main() {
    // Printing one computation
    if false {
        let mu = -0.01;
        let sigma = 0.1;
        let beta = 100.0;

        let alpha = -100.0;
        let parameters = Parameters::new(mu, sigma, alpha, beta).unwrap();
        let subs = Substitutions::new(N_SANDPIPER, parameters);

        println!("{:?}", subs.mean());
    }

    // Plotting the shape of fitting error while varying alpha
    if false {
        let mu = -0.004444444;
        let beta = 0.0;

        let sigma = 1.;
        let grid: Vec<f64> = ndarray::Array1::<f64>::linspace(-100., 100., 100).iter().map(|x| *x).collect();
        let values: Vec<f64> = grid.clone().iter()
            .map(|&alpha| {
                let parameters = Parameters::new(mu, sigma, alpha, beta).unwrap();
                let subs = Substitutions::new(N_SANDPIPER, parameters);
                (subs.mean() - EMPIRICAL_MEAN_SUBSTITUTIONS_SANDPIPER).abs().into()
            })
            .collect();

        let mut comparison: pre::Processes<_, _> = (grid.into_iter(), values)
            .preexplore()
            .set_title(format!("sigma = {}", sigma))
            .to_owned()
            .into();

        for sigma in [0.1, 0.01, 0.001, 0.0001].iter() {
            let grid: Vec<f64> = ndarray::Array1::<f64>::linspace(-100., 100., 100).iter().map(|x| *x).collect();
            let values: Vec<f64> = grid.clone().iter()
                .map(|&alpha| {
                    let parameters = Parameters::new(mu, *sigma, alpha, beta).unwrap();
                    let subs = Substitutions::new(N_SANDPIPER, parameters);
                    (subs.mean() - EMPIRICAL_MEAN_SUBSTITUTIONS_SANDPIPER).abs().into()
                })
                .collect();

            let line = (grid.into_iter(), values)
                .preexplore()
                .set_title(&format!("sigma = {}", sigma))
                .to_owned();

            comparison += line;
        }

        comparison
            .set_logy(2)
            .plot("error fitting for various alphas").unwrap();           

    }


    // Plotting the shape of fitting error while varying alpha, for various beta
    if false {
 
        let mu = -0.004444;
        let sigma = 0.1;
        let beta = 0.0;

        let grid: Vec<f64> = ndarray::Array1::<f64>::linspace(-300., 300., 100).iter().map(|x| *x).collect();
        let values: Vec<f64> = grid.clone().iter()
            .map(|&alpha| {
                let parameters = Parameters::new(mu, sigma, alpha, beta).unwrap();
                let subs = Substitutions::new(N_SANDPIPER, parameters);
                (subs.mean() - EMPIRICAL_MEAN_SUBSTITUTIONS_SANDPIPER).abs().into()
            })
            .collect();

        let mut comparison: pre::Processes<_,_> = (grid.into_iter(), values)
            .preexplore()
            .set_title(format!("beta = {}", beta))
            .set_xlabel("alpha")
            .set_ylabel("fitting error")
            .set_logy(2)
            .to_owned()
            .into();

        for beta in [1.0, 10.0, 100., 1000.].iter() {

            let grid: Vec<f64> = ndarray::Array1::<f64>::linspace(-300., 300., 100).iter().map(|x| *x).collect();
            let values: Vec<f64> = grid.clone().iter()
                .map(|&alpha| {
                    let parameters = Parameters::new(mu, sigma, alpha, *beta).unwrap();
                    let subs = Substitutions::new(N_SANDPIPER, parameters);
                    (subs.mean() - EMPIRICAL_MEAN_SUBSTITUTIONS_SANDPIPER).abs().into()
                })
                .collect();

            let line = (grid.into_iter(), values)
                .preexplore()
                .set_title(format!("beta = {}", beta))
                .set_xlabel("alpha")
                .set_ylabel("fitting error")
                .set_logy(2)
                .to_owned();

            comparison += line;
        }

        comparison
            .set_logy(2)
            .plot("alpha for various beta").unwrap();        

    
    }


    // Saving a heat map for each fixed beta.
    if true {
        for beta in [1000.].iter() { //, 10., 100., 200., 300., 500., 800., 1000., 1500., 2000.].iter() {
            fitting_inspection(*beta);
        }
    }
}

fn fitting_inspection(beta: f64) {
    // Params
    let grid_level = 10; 

    // Grids
    let mu_grid = ndarray::Array1::<f64>::linspace(-0.01, 0., grid_level);
    let sigma_grid = ndarray::Array1::<f64>::linspace(1. / 1000., 1. / 300., grid_level);

    // Searching for alpha
    let fitting_values: Vec<f64> = mu_grid
        .iter()
        .cartesian_product(sigma_grid.iter())
        .collect::<Vec<(&f64, &f64)>>()
        .into_par_iter()
        .map(|(mu, sigma)| {
            let alpha_fitting = |alpha: f64| -> f64 {
                let parameters = Parameters::new(*mu, *sigma, alpha, beta).unwrap();
                let subs = Substitutions::new(N_SANDPIPER, parameters);
                (subs.mean() - EMPIRICAL_MEAN_SUBSTITUTIONS_SANDPIPER)
                    .abs()
                    .into()
            };
            let derivative = |alpha: f64| -> f64 {
                let center = alpha_fitting(alpha);
                let right = alpha_fitting(alpha + std::f64::EPSILON);
                if center < right {
                    1.
                } else {
                    -1.
                }
            };
            let mut root_finding = ExpBinary::new(&derivative, 0.);
            root_finding.set_maximum(0.);
            match root_finding.run() {
                Ok(alpha) => alpha_fitting(alpha),
                Err(_) => std::f64::NAN,
            }
        })
        .collect();

    // Collect data
    let mut data = Vec::new();
    for i in 0..mu_grid.len() {
        for j in 0..sigma_grid.len() {
            data.push(mu_grid[i]);
            data.push(sigma_grid[j]);
            data.push(fitting_values[i * sigma_grid.len() + j].log(2.));
        }
    }

    // Plot it
    pre::Data::new(data, 3)
        .set_title(&format!("Fitting values for beta = {}", beta))
        .set_xlabel("mu")
        .set_ylabel("sigma")
        .plot_with_script(&format!("heat_map, beta = {}, sandpiper", beta), format!("
unset key
set title \"Fitting values for beta = {}\"
set xlabel \"mu\"
set ylabel \"sigma\"

set view map
plot \"target\\\\preexplorer\\\\data\\\\heat_map, beta = {}, sandpiper.txt\" u 1:2:3 with image
pause -1
            ", beta, beta)
        )
        .unwrap();
}