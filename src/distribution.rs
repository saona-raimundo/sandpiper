pub use self::beta::Beta;
pub use self::genetic_freq::GeneticFreq;
pub use self::heterozygosity::{Dominance, Heterozygosity, Selection};
pub use self::normal::Normal;
pub use self::skew_normal::SkewNormal;

mod beta;
mod genetic_freq;
mod heterozygosity;
mod normal;
mod skew_normal;
