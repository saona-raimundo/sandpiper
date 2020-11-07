// Traits
use rand::distributions::Distribution;
use rand::Rng;
use statrs::statistics::{Max, Min};

// Structs
use crate::error::{Result, StatsError};
use std::f64;

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Selection {
    Fixed(f64),
    SkewNormal{
        location: f64,
        scale: f64,
        shape: f64,
        bounds: Option<(f64, f64)>,
    }
}

#[derive(Debug, Copy, Clone, PartialEq)]
pub enum Dominance {
    /// Fixed value.
    Fixed(f64),
    Sigmoid{
        rate: f64,
    }
}

/// Polymorphisms in the population.
///
/// If 'x' is the allele frequency, then the heterozygosity is '2x(1-x)'.
///
/// # Examples
///
/// ```
/// use sandpiper::{Heterozygosity, Selection, Dominance};
///
/// let population = 1000;
/// let mutation_rate = 0.00001;
/// let selection = Selection::Fixed(0.0); 
/// let dominance = Dominance::Fixed(0.5);
/// let hetero = Heterozygosity::new(population, mutation_rate, selection, dominance).unwrap();
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Heterozygosity {
    population: u64, // N
    mutation_rate: f64, // U
    selection: Selection,
    dominance: Dominance,
}

impl Heterozygosity {
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
    /// use sandpiper::{Heterozygosity, Selection, Dominance};
    ///
    /// let population = 1000;
    /// let mutation_rate = 0.00001;
    /// let selection = Selection::Fixed(0.0); 
    /// let dominance = Dominance::Fixed(0.5);
    /// 
    /// let result = Heterozygosity::new(population, mutation_rate, selection, dominance);
    /// assert!(result.is_ok());
    ///
    /// let dominance = Dominance::Fixed(-0.5);
    /// let result = Heterozygosity::new(population, mutation_rate, selection, dominance);
    /// assert!(result.is_err());
    /// ```
    pub fn new(population: u64, mutation_rate: f64, selection: Selection, dominance: Dominance) -> Result<Self> {
        match selection {
            Selection::Fixed(s) => if s.is_nan() {
                return Err(StatsError::BadParams)
            },
            Selection::SkewNormal{location, scale, shape, bounds} => {
                if location.is_nan() || scale.is_nan() || scale <= 0.0 || shape.is_nan() {
                    return Err(StatsError::BadParams)
                }
                if let Some((lower_bound, upper_bound)) = bounds {
                    if upper_bound < lower_bound {
                        return Err(StatsError::BadParams)
                    }
                }
            },
        }
        match dominance {
            Dominance::Fixed(h) => if h.is_nan() || h < 0.0 {
                return Err(StatsError::BadParams)
            },
            Dominance::Sigmoid{rate} => if rate < 0.0 || rate.is_nan() {
                return Err(StatsError::BadParams)
            },
        }
        Ok(Heterozygosity{population, mutation_rate, selection, dominance})
    }

    /// Samples from the selection 's'.
    pub fn sample_selection<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        match self.selection {
            Selection::Fixed(s) => s,
            Selection::SkewNormal{location, scale, shape, bounds} => {
                let random_selection = crate::SkewNormal::new(location, scale, shape).unwrap();
                let mut proposal = random_selection.sample(rng);
                if let Some((lower_bound, upper_bound)) = bounds {
                    while (proposal < lower_bound) || (proposal > upper_bound) {
                        proposal = random_selection.sample(rng);
                    }
                }
                proposal
            },
        }
    }

    /// Samples from the allele frequency 'x' that will lead to heterozygosity '2 x (1 - x)'.
    pub fn sample_frequency<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let selection = self.sample_selection(rng);
        let dominance = match self.dominance {
            Dominance::Fixed(h) => h,
            Dominance::Sigmoid{rate} =>  1. / (1. + (-rate * selection).exp()),
        };

        crate::GeneticFreq::new(self.population, self.mutation_rate, selection, dominance)
            .unwrap()
            .sample(rng)
    }
}

impl Distribution<f64> for Heterozygosity {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> f64 {
        let x = self.sample_frequency(rng);

        2. * x * (1. - x)
    }
}

impl Min<f64> for Heterozygosity {
    /// Returns the minimum value in the domain
    fn min(&self) -> f64 {
        0.0
    }
}

impl Max<f64> for Heterozygosity {
    /// Returns the maximum value in the domain
    fn max(&self) -> f64 {
        0.5
    }
}

// #[cfg(test)]
// mod tests {
//     use rand::prelude::*;
//     use super::*;
// }