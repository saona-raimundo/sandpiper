// Traits
use rand::distributions::Distribution;
use rand::Rng;
use statrs::statistics::{Max, Min}; // , Mean, Variance};

// Structs
use crate::distribution::Beta;
use crate::error::{Result, StatsError};
use std::f64;

/// Distribution over allele frequency.
///
/// # Examples
///
/// ```
/// use sandpiper::GeneticFreq;
///
/// let population = 1000;
/// let mutation_rate = 0.00001;
/// let selection = -0.00001;
/// let dominance = 0.5;
/// let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance).unwrap();
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct GeneticFreq {
    population: u64,
    mutation_rate: f64,
    selection: f64,
    dominance: f64,
}

impl GeneticFreq {
    /// Constructs a new skew normal distribution with a location of `location`,
    /// a scale of `scale` and a shape of `shape`.
    ///
    /// # Errors
    ///
    /// Returns an error if `dominance` is not in the interval [0, 1]
    ///
    /// # Examples
    ///
    /// ```
    /// use sandpiper::GeneticFreq;
    ///
    /// let population = 1000;
    /// let mutation_rate = 0.00001;
    /// let selection = -0.00001;
    /// let dominance = 0.5;
    ///
    /// let result = GeneticFreq::new(population, mutation_rate, selection, dominance);
    /// assert!(result.is_ok());
    ///
    /// let dominance = -1.;
    /// let result = GeneticFreq::new(population, mutation_rate, selection, dominance);
    /// assert!(result.is_err());
    /// ```
    pub fn new(
        population: u64,
        mutation_rate: f64,
        selection: f64,
        dominance: f64,
    ) -> Result<Self> {
        if dominance.is_nan() || dominance > 1.0 || dominance < 0.0 {
            Err(StatsError::BadParams)
        } else {
            Ok(GeneticFreq {
                population,
                mutation_rate,
                selection,
                dominance,
            })
        }
    }
}

impl Distribution<f64> for GeneticFreq {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let shape = 4. * self.population as f64 * self.mutation_rate;
        let beta: Beta<f64> = Beta::new(shape, shape).unwrap();
        let mut proposal: f64 = beta.sample(rng);

        if self.selection < 0. {
            let reshaping = |x: f64| {
                (2. * self.population as f64
                    * self.selection
                    * (x.powi(2) + 2. * self.dominance * x * (1. - x)))
                    .exp()
            };
            while rng.sample::<f64, _>(rand_distr::Standard) > reshaping(proposal) {
                proposal = beta.sample(rng);
            }
        } else {
            let reshaping = |x: f64| {
                (2. * self.population as f64
                    * self.selection
                    * (x.powi(2) + 2. * self.dominance * x * (1. - x) - 1.))
                    .exp()
            };
            while rng.sample::<f64, _>(rand_distr::Standard) > reshaping(proposal) {
                proposal = beta.sample(rng);
            }
        }

        proposal
    }
}

impl Min<f64> for GeneticFreq {
    /// Returns the minimum value in the domain of the
    /// skew normal distribution representable by a double precision float
    ///
    /// # Formula
    ///
    /// ```ignore
    /// -INF
    /// ```
    fn min(&self) -> f64 {
        f64::NEG_INFINITY
    }
}

impl Max<f64> for GeneticFreq {
    /// Returns the maximum value in the domain of the
    /// skew normal distribution representable by a double precision float
    ///
    /// # Formula
    ///
    /// ```ignore
    /// INF
    /// ```
    fn max(&self) -> f64 {
        f64::INFINITY
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;

    #[test]
    fn mean() {
        let population = 1000;
        let mutation_rate = 0.00001;
        let selection = 0.0;
        let dominance = 0.5;

        let gen_freq = GeneticFreq::new(population, mutation_rate, selection, dominance).unwrap();
        let expected = 0.5;

        let sampled = gen_freq
            .sample_iter(crate::tests::rng(1))
            .take(100000)
            .sum::<f64>()
            / 100000.;

        assert!((expected - sampled).abs() < 1e-2);
    }
}
