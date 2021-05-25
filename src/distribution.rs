pub use self::beta::Beta;
pub use self::genetic_freq::GeneticFreq;
pub use self::heterozygosity::{
    Dominance, Heterozygosity, Selection, UnfixedHeterozygosity, UpperBound,
};
pub use self::normal::Normal;
pub use self::skew_normal::SkewNormal;

mod beta;
mod genetic_freq;
mod heterozygosity;
mod normal;
mod skew_normal;

pub mod helper {

    use rayon::prelude::*;

    /// Approximates the cummulative distribution of a random variable.
    ///
    /// # Algorithm
    ///
    /// Monte Carlo simulation until the error bound (between two samples) is met for the number of repetitions given.
    /// If this is not the case for the initial number of samples, the samples are doubled and freshly new sampled.
    /// The returned histogram contains all samples used along the computations
    ///
    /// # Remarks
    ///
    /// By the nature of the computation, the result is random. A proper error analysis can guarantee probabilistic bounds.
    ///
    /// # Examples
    ///
    /// ```
    /// ```
    pub fn approx_histogram<R: rand::Rng + ?Sized>(
        variable: impl rand_distr::Distribution<f64>,
        rng: &mut R,
        grid: Vec<f64>,
        init_samples: usize,
        repeteitions: usize,
        error: f64,
    ) -> anyhow::Result<quantiles::histogram::Histogram<f64>> {
        let mut samples = init_samples;
        let mut cum_histo = quantiles::histogram::Histogram::<f64>::new(grid.clone())
            .map_err(|_| anyhow::anyhow!("quantiles::histogram::Error"))?;
        let mut empirical_error = std::f64::INFINITY;
        loop {
            println!("Trying {} samples", samples);
            // First histogram
            // Simulation
            let mut init_histo = quantiles::histogram::Histogram::<f64>::new(grid.clone()).unwrap();
            for _ in 0..samples {
                init_histo.insert(variable.sample(rng));
            }
            // Recover
            cum_histo += init_histo.clone();
            // Distribution
            let init_distribution = (0..grid.len())
                .map(|i| {
                    init_histo.total_below(quantiles::histogram::Bound::Finite(grid[i])) as f64
                        / samples as f64
                })
                .collect::<Vec<f64>>();
            for _ in 0..repeteitions {
                // Other histogram
                let mut other_histo =
                    quantiles::histogram::Histogram::<f64>::new(grid.clone()).unwrap();
                for _ in 0..samples {
                    other_histo.insert(variable.sample(rng));
                }
                // Recover
                cum_histo += other_histo.clone();
                // Distribution
                let other_distribution = (0..grid.len())
                    .map(|i| {
                        other_histo.total_below(quantiles::histogram::Bound::Finite(grid[i])) as f64
                            / samples as f64
                    })
                    .collect::<Vec<f64>>();
                // Checking correctness
                empirical_error = (0..grid.len())
                    .map(|i| (init_distribution[i] - other_distribution[i]).abs())
                    .map(|x| ordered_float::NotNan::new(x).unwrap())
                    .max()
                    .unwrap()
                    .into_inner();
                if empirical_error > error {
                    break;
                }
            }
            if empirical_error < error {
                break;
            }
            samples *= 2;
        }
        Ok(cum_histo)
    }

    /// Approximates the cummulative distribution of a random variable with parallel sampling.
    ///
    /// # Algorithm
    ///
    /// Monte Carlo simulation until the error bound (between two samples) is met for the number of repetitions given.
    /// If this is not the case for the initial number of samples, the samples are doubled and freshly new sampled.
    /// The returned histogram contains only the last samples used.
    ///
    /// # Remarks
    ///
    /// By the nature of the computation, the result is random. A proper error analysis can guarantee probabilistic bounds.
    ///
    /// The rng used is `thread_rng`.
    ///
    /// # Examples
    ///
    /// ```
    /// ```
    pub fn par_approx_histogram(
        variable: impl rand_distr::Distribution<f64> + Sync,
        grid: Vec<f64>,
        init_samples: usize,
        repeteitions: usize,
        error: f64,
    ) -> anyhow::Result<quantiles::histogram::Histogram<f64>> {
        let mut samples = init_samples;
        loop {
            println!("Trying {} samples", samples);
            // First histogram
            // Simulation
            let mut init_histo = quantiles::histogram::Histogram::<f64>::new(grid.clone()).unwrap();
            let simulation: Vec<f64> = (0..samples)
                .into_par_iter()
                .map(|_| variable.sample(&mut rand::thread_rng()))
                .collect();
            for value in simulation {
                init_histo.insert(value);
            }
            // Distribution
            let init_distribution = (0..grid.len())
                .map(|i| {
                    init_histo.total_below(quantiles::histogram::Bound::Finite(grid[i])) as f64
                        / samples as f64
                })
                .collect::<Vec<f64>>();

            let result = (0..repeteitions).into_par_iter().find_any(|_| {
                // Other histogram
                let mut other_histo =
                    quantiles::histogram::Histogram::<f64>::new(grid.clone()).unwrap();
                let simulation: Vec<f64> = (0..samples)
                    .into_par_iter()
                    .map(|_| variable.sample(&mut rand::thread_rng()))
                    .collect();
                for value in simulation {
                    other_histo.insert(value);
                }
                // Distribution
                let other_distribution = (0..grid.len())
                    .map(|i| {
                        other_histo.total_below(quantiles::histogram::Bound::Finite(grid[i])) as f64
                            / samples as f64
                    })
                    .collect::<Vec<f64>>();
                // Checking correctness
                let empirical_error = (0..grid.len())
                    .map(|i| (init_distribution[i] - other_distribution[i]).abs())
                    .map(|x| ordered_float::NotNan::new(x).unwrap())
                    .max()
                    .unwrap()
                    .into_inner();
                // Checking
                empirical_error > error
            });
            match result {
                Some(_) => samples *= 2,
                None => break,
            }
        }
        // Final simulation
        let mut histo = quantiles::histogram::Histogram::<f64>::new(grid.clone()).unwrap();
        let simulation: Vec<f64> = (0..samples)
            .into_par_iter()
            .map(|_| variable.sample(&mut rand::thread_rng()))
            .collect();
        for value in simulation {
            histo.insert(value);
        }
        Ok(histo)
    }
}
