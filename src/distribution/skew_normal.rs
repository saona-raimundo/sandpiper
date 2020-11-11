// Traits
use rand::distributions::Distribution;
use rand::Rng;
use statrs::distribution::{Continuous, Univariate};
use statrs::statistics::{Max, Min}; // , Mean, Variance};

// Structs
use crate::error::{Result, StatsError};
use rand_distr::StandardNormal;
use std::f64;

/// Implements the [skew Normal](https://en.wikipedia.org/wiki/Skew_normal_distribution)
/// distribution.
///
/// # Examples
///
/// ```
/// use sandpiper::SkewNormal;
///
/// let sn = SkewNormal::new(0.0, 1.0, 0.0).unwrap();
/// ```
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct SkewNormal {
    location: f64,
    scale: f64,
    shape: f64,
}

impl SkewNormal {
    /// Constructs a new skew normal distribution with a location of `location`,
    /// a scale of `scale` and a shape of `shape`.
    ///
    /// # Errors
    ///
    /// Returns an error if `loc` or `scale` are `NaN` or if
    /// `scale <= 0.0`
    ///
    /// # Examples
    ///
    /// ```
    /// use sandpiper::SkewNormal;
    ///
    /// let mut result = SkewNormal::new(0.0, 1.0, 1.0);
    /// assert!(result.is_ok());
    ///
    /// result = SkewNormal::new(0.0, 0.0, 1.0);
    /// assert!(result.is_err());
    /// ```
    pub fn new(location: f64, scale: f64, shape: f64) -> Result<Self> {
        if location.is_nan() || scale.is_nan() || scale <= 0.0 {
            Err(StatsError::BadParams)
        } else {
            Ok(SkewNormal {
                location,
                scale,
                shape,
            })
        }
    }
}

impl Distribution<f64> for SkewNormal {
    fn sample<R: Rng + ?Sized>(&self, r: &mut R) -> f64 {
        sample_unchecked(r, self.location, self.scale, self.shape)
    }
}

impl Min<f64> for SkewNormal {
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

impl Max<f64> for SkewNormal {
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

impl Continuous<f64, f64> for SkewNormal {
    /// Calculates the probability density function for the
    /// skew normal distribution at `x`.
    ///
    /// # Formula
    ///
    /// ```ignore
    /// 2 / omega * phi( (x - xi) / omega ) * Phi( alpha (x - xi) / omega )
    /// ```
    ///
    /// where `xi` is the location, `omega` is the scale, `alpha` is the
    /// shape and `phi` and `Phi` are the density and distribution
    /// of a standard normal variable.
    fn pdf(&self, x: f64) -> f64 {
        pdf_unchecked(x, self.location, self.scale, self.shape)
    }

    /// Calculates the log probability density function for the skew normal
    /// distribution at `x`.
    ///
    /// # Formula
    ///
    /// ```ignore
    /// ln((1 / sqrt(2σ^2 * π)) * e^(-(x - μ)^2 / 2σ^2))
    /// ```
    ///
    /// where `μ` is the mean and `σ` is the standard deviation
    fn ln_pdf(&self, x: f64) -> f64 {
        ln_pdf_unchecked(x, self.location, self.scale, self.shape)
    }
}

/// performs an unchecked pdf calculation for a skew normal distribution
/// with the given mean and standard deviation at x
pub fn pdf_unchecked(x: f64, location: f64, scale: f64, shape: f64) -> f64 {
    let d = (x - location) / scale;
    let normal = crate::Normal::new(0., 1.).unwrap();
    2. / scale * normal.pdf(d) * normal.cdf(shape * d)
}

/// performs an unchecked log(pdf) calculation for a normal distribution
/// with the given mean and standard deviation at x
pub fn ln_pdf_unchecked(x: f64, location: f64, scale: f64, shape: f64) -> f64 {
    let d = (x - location) / scale;
    let normal = crate::Normal::new(0., 1.).unwrap();

    2.0_f64.ln() - scale.ln() + (-0.5 * d * d) - statrs::consts::LN_SQRT_2PI - scale.ln()
        + normal.cdf(shape * d).ln()
}

/// Draws two samples from a standard normal distribution using the
/// rand crate. Then, it follows the approach described by
/// [Ghorbanzadeh et. al.](http://dx.doi.org/10.4236/am.2014.513201)
pub fn sample_unchecked<R: Rng + ?Sized>(r: &mut R, location: f64, scale: f64, shape: f64) -> f64 {
    let linear_map = |x: f64| -> f64 { x * scale + location };
    let u_1: f64 = r.sample(StandardNormal);
    if shape == 0.0 {
        linear_map(u_1)
    } else {
        let u_2 = r.sample(StandardNormal);
        let (u, v) = (u_1.max(u_2), u_1.min(u_2));
        if shape == -1.0 {
            linear_map(v)
        } else if shape == 1.0 {
            linear_map(u)
        } else {
            let normalized = ((1. + shape) * u + (1. - shape) * v)
                / ((1. + shape * shape).sqrt() * std::f64::consts::SQRT_2);
            linear_map(normalized)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::prelude::*;
    use test_case::test_case;

    #[test_case(0., 1., 0., 0.0; "neutral shape")]
    #[test_case(0., 1., -1., -0.557; "negative shape")]
    #[test_case(0., 1., 1., 0.569; "positive shape")]
    #[test_case(0., 10., -1., -5.57; "spreaded negative shape")]
    #[test_case(0., 10., 1., 5.69; "spreaded positive shape")]
    fn mean(location: f64, scale: f64, shape: f64, expected: f64) {
        let skew_normal = SkewNormal::new(location, scale, shape).unwrap();
        let samples = 100_000;

        let result = skew_normal
            .sample_iter(crate::tests::rng(1))
            .take(samples)
            .sum::<f64>()
            / (samples as f64);

        println!("computed value: {:?}", result);
        assert!((expected - result).abs() < 1e-2);
    }
}
