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
mod distribution;
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
    #[test]
    fn it_works() {
        assert_eq!(2 + 2, 4);
    }
}
