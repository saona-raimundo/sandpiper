use rand::Rng;
use rand_distr::Float;
use rand_distr::{Distribution, Open01};
use std::{error, fmt};

/// The algorithm used for sampling the Beta distribution.
///
/// Reference:
///
/// R. C. H. Cheng (1978).
/// Generating beta variates with nonintegral shape parameters.
/// Communications of the ACM 21, 317-322.
/// https://doi.org/10.1145/359460.359482
#[derive(Clone, Copy, Debug)]
enum BetaAlgorithm<N> {
    BB(BB<N>),
    BC(BC<N>),
}

/// Algorithm BB for `min(alpha, beta) > 1`.
#[derive(Clone, Copy, Debug)]
struct BB<N> {
    alpha: N,
    beta: N,
    gamma: N,
}

/// Algorithm BC for `min(alpha, beta) <= 1`.
#[derive(Clone, Copy, Debug)]
struct BC<N> {
    alpha: N,
    beta: N,
    delta: N,
    kappa1: N,
    kappa2: N,
}

/// The Beta distribution with shape parameters `alpha` and `beta`.
///
/// # Example
///
/// ```
/// use rand_distr::{Distribution, Beta};
///
/// let beta = Beta::new(2.0, 5.0).unwrap();
/// let v = beta.sample(&mut rand::thread_rng());
/// println!("{} is from a Beta(2, 5) distribution", v);
/// ```
#[derive(Clone, Copy, Debug)]
pub struct Beta<N> {
    a: N,
    b: N,
    switched_params: bool,
    algorithm: BetaAlgorithm<N>,
}

/// Error type returned from `Beta::new`.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BetaError {
    /// `alpha <= 0` or `nan`.
    AlphaTooSmall,
    /// `beta <= 0` or `nan`.
    BetaTooSmall,
}

impl fmt::Display for BetaError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.write_str(match self {
            BetaError::AlphaTooSmall => "alpha is not positive in beta distribution",
            BetaError::BetaTooSmall => "beta is not positive in beta distribution",
        })
    }
}

impl error::Error for BetaError {}

impl<N: Float> Beta<N>
where
    Open01: Distribution<N>,
{
    /// Construct an object representing the `Beta(alpha, beta)`
    /// distribution.
    pub fn new(alpha: N, beta: N) -> Result<Beta<N>, BetaError> {
        if !(alpha > N::from(0.)) {
            return Err(BetaError::AlphaTooSmall);
        }
        if !(beta > N::from(0.)) {
            return Err(BetaError::BetaTooSmall);
        }
        // From now on, we use the notation from the reference,
        // i.e. `alpha` and `beta` are renamed to `a0` and `b0`.
        let (a0, b0) = (alpha, beta);
        let (a, b, switched_params) = if a0 < b0 {
            (a0, b0, false)
        } else {
            (b0, a0, true)
        };
        if alpha > N::from(1.) {
            let alpha = a + b;
            let beta = ((alpha - N::from(2.)) / (N::from(2.) * a * b - alpha)).sqrt();
            let gamma = a + N::from(1.) / beta;

            Ok(Beta {
                a,
                b,
                switched_params,
                algorithm: BetaAlgorithm::BB(BB { alpha, beta, gamma }),
            })
        } else {
            let alpha = a + b;
            let beta = N::from(1.) / b;
            let delta = N::from(1.) + a - b;
            let kappa1 = delta * (N::from(1. / 18. / 4.) + N::from(3. / 18. / 4.) * b)
                / (a * beta - N::from(14. / 18.));
            let kappa2 = N::from(0.25) + (N::from(0.5) + N::from(0.25) / delta) * b;

            Ok(Beta {
                a,
                b,
                switched_params,
                algorithm: BetaAlgorithm::BC(BC {
                    alpha,
                    beta,
                    delta,
                    kappa1,
                    kappa2,
                }),
            })
        }
    }
}

impl<N: Float> Distribution<N> for Beta<N>
where
    Open01: Distribution<N>,
{
    fn sample<R: Rng + ?Sized>(&self, rng: &mut R) -> N {
        match self.algorithm {
            BetaAlgorithm::BB(algo) => {
                let mut w;
                loop {
                    // 1.
                    let u1 = rng.sample(Open01);
                    let u2 = rng.sample(Open01);
                    let v = algo.beta * (u1 / (N::from(1.) - u1)).ln();
                    w = self.a * v.exp();
                    let z = u1 * u1 * u2;
                    let r = algo.gamma * v - N::from(4.).ln();
                    let s = self.a + r - w;
                    // 2.
                    if s + N::from(1.) + N::from(5.).ln() >= N::from(5.) * z {
                        break;
                    }
                    // 3.
                    let t = z.ln();
                    if s >= t {
                        break;
                    }
                    // 4.
                    if !(r + algo.alpha * (algo.alpha / (self.b + w)).ln() < t) {
                        break;
                    }
                }
                // 5.
                if !self.switched_params {
                    w / (self.b + w)
                } else {
                    self.b / (self.b + w)
                }
            }
            BetaAlgorithm::BC(algo) => {
                let mut w;
                loop {
                    let z;
                    // 1.
                    let u1 = rng.sample(Open01);
                    let u2 = rng.sample(Open01);
                    if u1 < N::from(0.5) {
                        // 2.
                        let y = u1 * u2;
                        z = u1 * y;
                        if N::from(0.25) * u2 + z - y >= algo.kappa1 {
                            continue;
                        }
                    } else {
                        // 3.
                        z = u1 * u1 * u2;
                        if z <= N::from(0.25) {
                            let v = algo.beta * (u1 / (N::from(1.) - u1)).ln();
                            w = self.a * v.exp();
                            break;
                        }
                        // 4.
                        if z >= algo.kappa2 {
                            continue;
                        }
                    }
                    // 5.
                    let v = algo.beta * (u1 / (N::from(1.) - u1)).ln();
                    w = self.a * v.exp();
                    if !(algo.alpha * ((algo.alpha / (self.b + w)).ln() + v) - N::from(4.).ln()
                        < z.ln())
                    {
                        break;
                    };
                }
                // 6.
                if !self.switched_params {
                    if w == N::from(std::f64::INFINITY) {
                        // Assuming `b` is finite, for large `w`:
                        return N::from(1.);
                    }
                    w / (self.b + w)
                } else {
                    self.b / (self.b + w)
                }
            }
        }
    }
}
