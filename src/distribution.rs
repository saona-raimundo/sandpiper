pub use self::beta::Beta;
pub use self::genetic_freq::GeneticFreq;
pub use self::normal::Normal;
pub use self::skew_normal::SkewNormal;
pub use self::heterozygosity::{Heterozygosity, Selection, Dominance};

mod genetic_freq;
mod normal;
mod skew_normal;
mod beta;
mod heterozygosity;