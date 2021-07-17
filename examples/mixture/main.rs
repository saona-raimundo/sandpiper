use itertools::iproduct;
use preexplorer::prelude::*;
use read_input::prelude::*;
use sandpiper::prelude::*;
use splines::Spline;
use statrs::distribution::Continuous;

const BETA: f64 = 0.;
const MU_MIN: f64 = -0.1;
const MU_POINTS: i32 = 10;
const POPULATION_SIZE: u64 = 500_000;
const UPPER_GEN_FREQ: sandpiper::UpperBound = sandpiper::UpperBound::Smallest;
const VARIANCE_SAMPLES: usize = 1_000;
const ERROR_LIMIT: f64 = 1e-3;
const SIGMA_MAX: f64 = 0.1;
const SIGMA_POINTS: i32 = 10;

fn main() -> anyhow::Result<()> {
    // Initialization
    println!(
        "\
Hello! Welcome!
This program will compute interatively \
expected polymorphism for fixed beta and the case of \
alpha = 0 (ie no skewness, just a normal gaussian distribution)."
    );
    println!("Let's confirm some values.");

    let beta: f64 = input()
        .msg(format!("Please input beta (default {}): ", BETA))
        .default(BETA)
        .get();
    let population_size = input()
        .msg(format!(
            "Please input population size (default {}): ",
            POPULATION_SIZE
        ))
        .default(POPULATION_SIZE)
        .get();
    let mu_min = input()
        .msg(format!("Please input minimum mu (default {}): ", MU_MIN))
        .default(MU_MIN)
        .add_err_test(|&x| x < 0., "Please input a negative value.")
        .add_err_test(|&x| x > -1., "Please input a value greater than -1.")
        .get();
    let mu_points = input()
        .msg(format!(
            "Please input number of points in the grid (default {}): ",
            MU_POINTS
        ))
        .default(MU_POINTS)
        .add_err_test(|&x| x > 0, "Please input a positive value.")
        .get();

    // Computing initial conditions (fixed selection)
    let initial_conditions: Spline<f64, f64>;
    {
        println!("Computing a symetric geometric grid to zero.");
        let grid_geometric = ndarray::Array::geomspace(
            mu_min / (2.0_f64.powi(mu_points)),
            mu_min,
            mu_points as usize,
        )
        .ok_or(anyhow::anyhow!("geometric grid failed!"))?;

        let mut grid_initial = grid_geometric.clone().to_vec();
        grid_initial.push(0.);
        grid_initial.extend(grid_geometric.clone().into_iter().map(|mu| mu.abs()));

        let conditions_initial = grid_initial.iter().map(|&location| {
            let hetero = UnfixedHeterozygosity::new(
                population_size,
                sandpiper::U,
                Selection::Fixed(location),
                Dominance::Sigmoid { rate: beta },
                UPPER_GEN_FREQ,
            )
            .unwrap();
            hetero.mc_approx_mean(VARIANCE_SAMPLES, ERROR_LIMIT)
        });

        initial_conditions = Spline::from_iter(grid_initial.iter().zip(conditions_initial).map(
            |(&location, value)| {
                splines::Key::new(location, value.mean(), splines::Interpolation::Linear)
            },
        ));

        // Plotting initial conditions
        if true {
            initial_conditions
                .keys()
                .iter()
                .map(|key| (key.t, key.value))
                .unzip::<_, _, Vec<f64>, Vec<f64>>()
                .preexplore()
                .set_ylog(10)
                .plot("grid")?;
        }
    }

    let sigma_max = input()
        .msg(format!(
            "Please input maximum sigma (default {}): ",
            SIGMA_MAX
        ))
        .default(SIGMA_MAX)
        .add_err_test(|&x| x > 0., "Please input a negative value.")
        .get();
    let sigma_points = input()
        .msg(format!(
            "Please input number of points (default {}): ",
            SIGMA_POINTS
        ))
        .default(SIGMA_POINTS)
        .add_err_test(|&x| x > 0, "Please input a positive value.")
        .get();
    // Iteration
    let mut temporal_values = vec![initial_conditions.clone()];
    let grid_geometric = ndarray::Array::geomspace(
        sigma_max / (2.0_f64.powi(sigma_points)),
        sigma_max,
        sigma_points as usize,
    )
    .ok_or(anyhow::anyhow!("geometric grid failed!"))?;
    {
        println!("Computing iterations!");

        for temporal_value in grid_geometric.iter() {
            let mut new_keys = vec![];

            // Convolution for each point
            for spatial_index in 0..initial_conditions.len() {
                let spatial_center: f64 = initial_conditions
                    .get(spatial_index)
                    .ok_or(anyhow::anyhow!(format!(
                        "There is no key #{}!",
                        spatial_index
                    )))
                    .unwrap()
                    .t;
                // Convolution
                let gaussian = statrs::distribution::Normal::new(0.0, *temporal_value).unwrap();
                let integrand = |x: f64| -> f64 {
                    gaussian.pdf(spatial_center - x) * initial_conditions.clamped_sample(x).unwrap()
                };
                let new_value = quadrature::integrate(
                    integrand,
                    spatial_center - 4. * temporal_value,
                    spatial_center + 4. * temporal_value,
                    1e-10,
                )
                .integral;
                // Retrieving the value
                new_keys.push(splines::Key::new(
                    spatial_center,
                    new_value,
                    splines::Interpolation::Linear,
                ));
            }
            temporal_values.push(Spline::from_vec(new_keys));
        }
    }

    // Plotting result
    let (grid, _) = initial_conditions
        .keys()
        .iter()
        .map(|key| (key.t, key.value))
        .unzip::<_, _, Vec<f64>, Vec<f64>>();

    let values =
        iproduct!(grid.iter(), 0..SIGMA_POINTS as usize).map(|(spatial_value, temporal_index)| {
            temporal_values[temporal_index]
                .clamped_sample(*spatial_value)
                .unwrap()
        });
    pre::Heatmap::new(grid.iter(), grid_geometric.iter(), values)
        .plot("simple")
        .unwrap();

    Ok(())
}
