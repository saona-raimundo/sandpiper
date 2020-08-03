use crate::distribution::SkewNormal;
use crate::error::Result;

/// Parameters of the model.
///
#[derive(Debug, Copy, Clone)]
pub struct Parameters {
    pub mu: f64,
    pub sigma: f64,
    pub alpha: f64,
    pub beta: f64,
    pub skew_normal: SkewNormal,
}

impl Parameters {
    pub fn new(mu: f64, sigma: f64, alpha: f64, beta: f64) -> Result<Self> {
        let skew_normal = SkewNormal::new(mu, sigma, alpha)?;
        Ok(Parameters {
            mu,
            sigma,
            alpha,
            beta,
            skew_normal,
        })
    }
}
