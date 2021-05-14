//! Sandpiper computations helper trait.
//!
//! This crates includes helper functions for computations used in the sandpiper project.

pub use self::constants::*;
pub use self::distribution::*;
pub use self::parameters::Parameters;
pub use self::root_finding::ExpBinary;
pub use self::statistic::*;

/// Empirical data and overall constants.
mod constants;
/// Distributions.
pub mod distribution;
/// Errors and results from this crate.
pub mod error;
/// Parameters of the model.
mod parameters;
/// Root finding algorithms.
mod root_finding;
/// Statistics of concern in the sandpiper.
mod statistic;

pub mod prelude {
    pub use crate::constants::*;
    pub use crate::distribution::*;
    pub use crate::{Parameters, Substitutions};
    pub use statrs::statistics::Mean;
}

#[cfg(test)]
mod tests {
    /// Construct a deterministic RNG with the given seed
    pub fn rng(seed: u64) -> impl rand::RngCore {
        // For tests, we want a statistically good, fast, reproducible RNG.
        // PCG32 will do fine, and will be easy to embed if we ever need to.
        const INC: u64 = 11634580027462260723;
        rand_pcg::Pcg32::new(seed, INC)
    }

    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
