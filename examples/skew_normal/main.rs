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
        let locations = vec![
            -0.0000600000,
            -0.0000566291,
            -0.0000534476,
            -0.0000504448,
            -0.0000476107,
            -0.0000449359,
            -0.0000424113,
            -0.0000400285,
            -0.0000377797,
            -0.0000356571,
            -0.0000336539,
            -0.0000317631,
            -0.0000299786,
            -0.0000282944,
            -0.0000267047,
            -0.0000252044,
            -0.0000237884,
            -0.0000224519,
            -0.0000211905,
            -0.0000200000,
        ]; // vec![-4e-5]; //
        let scales = vec![
            0.0000100000,
            0.0000117078,
            0.0000137073,
            0.0000160482,
            0.0000187889,
            0.0000219977,
            0.0000257544,
            0.0000301527,
            0.0000353022,
            0.0000413311,
            0.0000483897,
            0.0000566536,
            0.0000663290,
            0.0000776566,
            0.0000909188,
            0.0001064459,
            0.0001246247,
            0.0001459081,
            0.0001708263,
            0.0002000000,
        ]; // vec![1e-5]; //
        let shapes = vec![-5.];

        let mut densities_vec = Vec::new();
        for location in locations {
            for scale in &scales {
                for shape in &shapes {
                    let sn = SkewNormal::new(location, *scale, *shape)?;
                    let realizations = (0..SAMPLES).map(|_| sn.sample(&mut rand::thread_rng()));
                    densities_vec.push(
                        pre::Density::new(realizations)
                            .set_title(format!("({}, {}, {})", location, scale, shape))
                            .to_owned(),
                    )
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
