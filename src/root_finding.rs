use anyhow::{anyhow, Result};
use getset::{Getters, Setters};

/// Root finding method with exponential search for space together with binary search for root.
///
/// Inputs `maximum` and `minimum` bound the search space. They could be infinity, meaning that
/// there is no known bound on the space in which the root might be. From an intial point,
/// the algorithm expands the lower and upper bounds in both directions at exponential rate.
/// Once an interval where the lower bound evaluates a different sign from the upper bound
/// is found, the binary search starts.
#[derive(Debug, Getters, Setters)]
pub struct ExpBinary<F> {
    /// Function to find the root of.
    #[getset(set = "pub", get = "pub")]
    f: F,

    /// Tolerance when searching the root.
    ///
    /// This tolerance refers to the size of the interval in which the root is to be founded.
    #[getset(set = "pub", get = "pub")]
    tol: f64,

    /// Overall maximum bound.
    #[getset(set = "pub", get = "pub")]
    maximum: f64,

    /// Overall minimum bound.
    #[getset(set = "pub", get = "pub")]
    minimum: f64,

    /// Maximum number of steps in the exponential search.
    #[getset(set = "pub", get = "pub")]
    max_exp_iterations: u64,

    /// Step for exponential search.
    ///
    /// It is multiplied by 2 at each iteration.
    #[getset(set = "pub", get = "pub")]
    step: f64,

    counter_exp_iteration: u64,
    upper: f64,
    lower: f64,
}

impl<F> ExpBinary<F>
where
    F: Fn(f64) -> f64,
{
    /// Constructor
    pub fn new(f: F, init: f64) -> Self {
        // Dafualt values
        let tol = 1e-10;
        let maximum = std::f64::INFINITY;
        let minimum = -std::f64::INFINITY;
        let max_exp_iterations = 1000;
        let counter_exp_iteration = 0;
        let step = 1e-6;

        let upper = init;
        let lower = init;
        ExpBinary {
            f,
            tol,
            minimum,
            maximum,
            max_exp_iterations,
            counter_exp_iteration,
            step,
            upper,
            lower,
        }
    }

    pub fn run(&mut self) -> Result<f64> {
        // Exponential search
        self.upper = (self.upper + self.step).min(self.maximum);
        self.lower = (self.lower - self.step).max(self.minimum);
        self.step *= 2.;

        loop {
            // First step verification
            if (self.f)(self.lower) * (self.f)(self.upper) < 0. {
                break;
            }
            // Upper search
            let next_upper = (self.upper + self.step).min(self.maximum);
            if (self.f)(self.upper) * (self.f)(next_upper) < 0. {
                self.lower = self.upper;
                self.upper = next_upper;
                break;
            }
            self.upper = next_upper;

            // Lower search
            let next_lower = (self.lower - self.step).max(self.minimum);
            if (self.f)(self.lower) * (self.f)(next_lower) < 0. {
                self.upper = self.lower;
                self.lower = next_lower;
                break;
            }
            self.lower = next_lower;

            // Update
            self.step *= 2.;
            self.counter_exp_iteration += 1;

            if self.counter_exp_iteration > self.max_exp_iterations {
                return Err(anyhow!(
                    "Early stopped exponential search: {} iterations reached",
                    self.counter_exp_iteration
                ));
            }
        }
        assert!((self.f)(self.upper) * (self.f)(self.lower) <= 0.0);

        // Binary search
        let mut upper_value = (self.f)(self.upper);

        while (self.lower - self.upper).abs() > self.tol {
            let middle = (self.lower + self.upper) * 0.5;
            let middle_value = (self.f)(middle);

            if middle_value * upper_value > 0. {
                self.upper = middle;
                upper_value = middle_value;
            } else {
                self.lower = middle;
            }
        }
        Ok((self.lower + self.upper) * 0.5)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use test_case::test_case;

    #[test_case(|x| x, 0. ; "identity")]
    #[test_case(|x| x, 100. ; "identity from positive far away")]
    #[test_case(|x| x, -100. ; "identity from negative far away")]
    #[test_case(|x| x * (x - 100.), -100. ; "quadratic from negative far away")]
    #[test_case(|x| x * (x - 1000.), 200. ; "quadratic from positive far away")]
    fn trivial_functions<F>(f: F, init: f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut root_finding = ExpBinary::new(&f, init);
        let root = root_finding.run().unwrap();
        println!("Root finded: {}", root);
        assert!(root.abs() < 1e-6);
    }

    #[test_case(|x| x, 0. ; "identity")]
    #[test_case(|x| x, 100. ; "identity from positive far away")]
    #[test_case(|x| x, -100. ; "identity from negative far away")]
    #[test_case(|x| x * (x - 100.), -100. ; "quadratic from negative far away")]
    #[test_case(|x| x * (x - 1000.), 200. ; "quadratic from positive far away")]
    fn tolerance_control<F>(f: F, init: f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut root_finding = ExpBinary::new(&f, init);
        root_finding.set_tol(1e-15);
        println!("{:?}", root_finding.tol());

        let root = root_finding.run().unwrap();
        println!("Root finded: {}", root);
        assert!(root.abs() < 1e-14);
    }

    #[test_case(|x| x * (x - 1000.), 900. ; "quadratic from positive far away")]
    fn maximum_control<F>(f: F, init: f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut root_finding = ExpBinary::new(&f, init);
        root_finding.set_maximum(900.);
        println!("{:?}", root_finding.maximum());

        let root = root_finding.run().unwrap();
        println!("Root finded: {}", root);
        assert!(root < 1e-6);
    }

    #[test_case(|x| { if x < 0. { -(-x).sqrt() } else { x.sqrt() } }, 0. ; "symmetric sqrt")]
    #[test_case(|x| { if x < 0. { -(-x).sqrt() } else { x.sqrt() } }, 100. ; "symmetric sqrt from positive far away")]
    #[test_case(|x| { if x < 0. { -(-x).sqrt() } else { x.sqrt() } }, -100. ; "symmetric sqrt from negative far away")]
    fn sharp_functions<F>(f: F, init: f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut root_finding = ExpBinary::new(&f, init);
        let root = root_finding.run().unwrap();
        println!("Root finded: {}", root);
        assert!(root.abs() < 1e-6);
    }

    #[test_case(|x| { if x < 0. { -(-x).sqrt() - 1. } else { x.sqrt() + 1. } }, 0. ; "displaced symmetric sqrt")]
    #[test_case(|x| { if x < 0. { -(-x).sqrt() - 1. } else { x.sqrt() + 1. } }, 100. ; "displaced symmetric sqrt from positive far away")]
    #[test_case(|x| { if x < 0. { -(-x).sqrt() - 1. } else { x.sqrt() + 1. } }, -100. ; "displaced symmetric sqrt from negative far away")]
    fn non_root_functions<F>(f: F, init: f64)
    where
        F: Fn(f64) -> f64,
    {
        let mut root_finding = ExpBinary::new(&f, init);
        let root = root_finding.run().unwrap();
        println!("Root finded: {}", root);
        assert!(root.abs() < 1e-6);
    }
}
