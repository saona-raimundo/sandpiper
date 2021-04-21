mod constants {
    // Model parameters
    pub const MUS: [f64; 1] = [-0.0195556029];
    pub const SIGMAS: [f64; 1] = [0.0132709376];
    pub const ALPHAS: [f64; 1] = [0.];
    pub const BETAS: [f64; 1] = [0.];
    pub const TOTAL: usize = MUS.len() * SIGMAS.len() * ALPHAS.len() * BETAS.len();

    // Simulation parameters
    pub const LOWER_S: f64 = -1.;
    pub const UPPER_S: f64 = 1.;
    pub const VARIANCE_SAMPLES: usize = 1_000;
    pub const ERROR_LIMIT: f64 = 1e-6;
}

mod collective;
mod individual;


fn main() {
    if true {
        collective::collective_main();
    } else {
        individual::individual_main();
    }
}