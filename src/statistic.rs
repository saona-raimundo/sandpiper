//! Statistics that can be contrasted with real data.

// Crates
use noisy_float::prelude::*;
use quadrature::integrate;

// Traits
use rand::distributions::Distribution;
use rand::Rng;
use statrs::{
    distribution::{Continuous, Univariate},
    statistics::Mean,
};

// Structs
use crate::distribution::Normal;
use crate::Parameters;

// Types
use noisy_float::types::R64; // finite f64.

// Constants
use crate::constants::{EPS, T, U};

#[derive(Debug, Copy, Clone)]
pub struct Substitutions {
    n: R64,
    u: R64,
    t: R64,
    param: Parameters,
}

impl Substitutions {
    pub fn new(n: u64, param: Parameters) -> Self {
        let n: R64 = r64(n as f64);
        let t: R64 = r64(T as f64);
        let u: R64 = r64(U);

        Substitutions { n, u, t, param }
    }

    fn h(&self, s: R64) -> R64 {
        if s >= 0.0 {
            r64(0.5)
        } else {
            r64(0.5) * (s * self.param.beta).exp()
        }
    }

    /// Returns the number of mutations in one generation and in case of mutation
    /// for the given positive selection coefficient.
    ///
    fn positive_selection(&self, n: R64, s: R64) -> R64 {
        let partial;
        if s < 1e-6 {
            let n_inv = r64(1.0) / n;
            partial = n_inv * 0.5
                + (-n_inv * 0.5 + 1.) * 0.5 * s
                + (n * 2. + n_inv - 3.) * 0.08333333333333333 * s * s
        } else {
            partial = (-(-s).exp() + 1.) / (-(-n * 2. * s).exp() + 1.)
        }
        n * partial
    }

    /// Returns the number of mutations in one generation and in case of mutation
    /// for the given negative selection coefficient.
    ///
    fn negative_selection(&self, n: R64, s: R64) -> f64 {
        let h = self.h(s);

        let numerator_integrand = |x: f64| -> f64 {
            (-s * 2. * ((-h * 2. + 1.) / n * x * x + h * 2. * x - (-h * 2. + 1.) / (n * 4.) - h))
                .exp()
                .into()
        };
        let denominator_integrand = |x: f64| -> f64 {
            (-n * 2. * s * ((-h * 2. + 1.) * x * x + h * 2. * x - 1.))
                .exp()
                .into()
        };
        let scale: f64 = (-s * 2. * ((-h * 2. + 1.) / 4. * n + h - n)).exp().into();

        let numerator = integrate(numerator_integrand, 0.0, 0.5, EPS).integral;
        let denominator = integrate(denominator_integrand, 0.0, 1.0, EPS).integral;

        numerator / denominator * scale
    }
}

impl Distribution<R64> for Substitutions {
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> R64 {
        let s = r64(self.param.skew_normal.sample(rng));
        let h = self.h(s);

        let aux = h * -2. + 1.;

        if aux == 0. {
            ((-s).exp() - 1.) / ((self.n * s * -2.).exp() - 1.)
        } else {
            let normal = Normal::new(0., 1.).unwrap();

            let current_mean = -h / aux;

            println!("{:?}", self.n);
            println!("{:?}", s);
            println!("{:?}", (h * -2. + 1.));
            println!("{:?}", self.n * s * (h * -2. + 1.));
            let current_inv_std = (self.n * s * (h * -2. + 1.)).sqrt() * 2.;

            let normalize = |x: R64| -> R64 { (x - current_mean) * current_inv_std };
            let value = normal.cdf(normalize(r64(1.) / (self.n * 2.)).into())
                - normal.cdf(normalize(r64(0.)).into());

            self.n * self.u * self.t * value * 2.
                / (normal.cdf(normalize(r64(1.)).into()) - normal.cdf(normalize(r64(0.)).into()))
        }
    }
}

impl Mean<R64> for Substitutions {
    /// Returns the expected value of substitutions.
    ///
    /// # Examples
    ///
    /// ```
    /// use sandpiper::prelude::*;
    ///
    /// let mu = -0.01;
    /// let sigma = 1.0;
    /// let alpha = -1.0;
    /// let beta = 1000.0;
    /// let parameters = Parameters::new(mu, sigma, alpha, beta).unwrap();
    /// let subs = Substitutions::new(N_SANDPIPER, parameters);
    ///
    /// println!("Expected value of substitutions: {:?}", subs.mean());
    /// ```
    fn mean(&self) -> R64 {
        let integrand = |s: f64| -> f64 {
            if s > 0.0 {
                (self.positive_selection(self.n, r64(s)) * self.param.skew_normal.pdf(s)).into()
            } else {
                (self.negative_selection(self.n, r64(s)) * self.param.skew_normal.pdf(s)).into()
            }
        };
        let expectation = integrate(integrand, -10.0, 10.0, EPS).integral;

        self.u * self.t * expectation * 2.
    }
}


#[cfg(test)]
mod tests {
    use super::*;
    use test_case::test_case;

    #[test_case(0., 100., 0.5; "zero")]
    #[test_case(1., 100., 0.5; "positive")]
    #[test_case(10., 100., 0.5; "large positive")]
    #[test_case(-1., 100., 0.5 * (-100.).exp(); "negative")]
    #[test_case(-10., 100., 0.5 * (-1000.).exp(); "large negative")]
    fn dominance_coefficient_h(s: f64, beta: f64, expected: f64) {
        for mu in [-1., 0., 1.].iter() {
            for sigma in [0.1, 1., 2.].iter() {
                for alpha in [-1., 0., 1.].iter() {
                    let parameters = Parameters::new(*mu, *sigma, *alpha, beta).unwrap();
                    let subs = Substitutions::new(crate::N_SANDPIPER, parameters);

                    assert_eq!(subs.h(r64(s)), expected);        
                }
            }
        }
    }
}