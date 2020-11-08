use sandpiper::prelude::*;

fn main() {
	let hetero = Heterozygosity::new(
		N_REDNECK,
		U,
		Selection::Fixed(0.0),
		Dominance::Fixed(0.5),
		)
		.unwrap();

	println!("{:?}", hetero.mc_approx_mean(1000, 1e-5));


	// let hetero = Heterozygosity::new(
	// 	500000,
	// 	U,
	// 	Selection::SkewNormal{location: -0.00001, scale: 0.0001, shape: -2., bounds: None},
	// 	Dominance::Fixed(0.5),
	// 	)
	// 	.unwrap();

	// println!("{:?}", hetero.mc_approx_mean(1000, 1e-6));



}