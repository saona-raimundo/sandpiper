use anyhow::Result;
use preexplorer::prelude::*;
use rand::distributions::Distribution;
use sandpiper::prelude::*;

const SAMPLES: u64 = 10_000;

fn main() -> Result<()> {
    // Single skew normal
    if false {
        let location = -1.;
        let scale = 1.;
        let shape = 2.;
        let sn = SkewNormal::new(location, scale, shape)?;
        let mut rng = rand::thread_rng();

        let realizations = (0..SAMPLES).map(|_| sn.sample(&mut rng));

        pre::Density::new(realizations)
            .set_title("Skew normal")
            .plot("testing")
            .unwrap();
    }

    // Single Normal
    if false {
        let n = Normal::new(-1., 2.)?;
        let mut rng = rand::thread_rng();

        let realizations = (0..SAMPLES).map(|_| n.sample(&mut rng));

        pre::Density::new(realizations)
            .set_title("Normal")
            .plot("testing")
            .unwrap();
    }

    // Various Skew normals
    if true {
        let locations = vec![-6e-5, -4e-5, -2e-5];
        let scales = vec![1e-5, 2e-5];
        let shapes = vec![-8., -5., -2.];

        let mut densities_vec = Vec::new();
        for location in locations {
            for scale in &scales {
                for shape in &shapes {
                    let sn = SkewNormal::new(location, *scale, *shape)?;
                    let realizations = (0..SAMPLES).map(|_| sn.sample(&mut rand::thread_rng()));
                    densities_vec.push(pre::Density::new(realizations).set_title(format!("({}, {}, {})", location, scale, shape)).to_owned())
                }
            }
        }

        pre::Densities::new(densities_vec)
            .set_title("Skew normals of interest")
            .set_xlabel("selection coefficient")
            .set_ylabel("value")
            .plot("skew_normals")
            .unwrap();
    }

    Ok(())
}
