use itertools::iproduct;
use preexplorer::prelude::*;
use read_input::prelude::*;
use sandpiper::prelude::*;
use splines::Spline;

const BETA: f64 = 0.;
const MU_MIN: f64 = -0.1;
const MU_POINTS: i32 = 10;
const POPULATION_SIZE: u64 = 500_000;
const UPPER_GEN_FREQ: sandpiper::UpperBound = sandpiper::UpperBound::Smallest;
const VARIANCE_SAMPLES: usize = 1_000;
const ERROR_LIMIT: f64 = 1e-3;
const ITERATIONS: usize = 1000;

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
        .msg(format!("Please input minimum mu (default {}): ", MU_POINTS))
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

    let dx = mu_min.abs() / (2.0_f64.powi(mu_points)); // mu_min.abs() / mu_points.pow(2) as f64;
    let grid = ndarray::Array::range(mu_min, mu_min.abs(), dx);
    let mut temporal_values: Vec<Vec<f64>> = vec![grid
        .iter()
        .map(|mu| {
            initial_conditions
                .sample(*mu)
                .ok_or(anyhow::anyhow!(format!("Could not sample {}!", mu)))
                .unwrap()
        })
        .collect()];
    // Iteration
    let iterations = input()
        .msg(format!(
            "Please input number of iterations (default {}): ",
            ITERATIONS
        ))
        .default(ITERATIONS)
        .get();
    for temporal_index in 0..iterations {
        let mut last = vec![0.; grid.len()];
        last[0] =
            (3. * temporal_values[temporal_index][0] + temporal_values[temporal_index][1]) / 4.;
        for spatial_index in 1..(grid.len() - 1) {
            last[spatial_index] = (temporal_values[temporal_index][spatial_index - 1]
                + 2. * temporal_values[temporal_index][spatial_index]
                + temporal_values[temporal_index][spatial_index + 1])
                / 4.;
        }
        last[grid.len() - 1] = (temporal_values[temporal_index][grid.len() - 2]
            + 3. * temporal_values[temporal_index][grid.len() - 1])
            / 4.;
        temporal_values.push(last);
    }

    // Plotting result
    let values = iproduct!(0..grid.len(), 0..iterations)
        .map(|(spatial_index, temporal_index)| temporal_values[temporal_index][spatial_index]);
    pre::Heatmap::new(grid.iter(), 0..iterations, values)
        .plot("simple")
        .unwrap();

    Ok(())
}
