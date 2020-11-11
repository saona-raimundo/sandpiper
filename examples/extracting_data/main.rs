use std::fs::File;
use std::io::{BufRead, BufReader, Error, Write};

fn main() -> Result<(), Error> {
    // let target_path = "target/preexplorer/data/sandpiper_poly_all.txt";
    let target_path = "target/preexplorer/data/redneck_poly_all.txt";
    let mut output = File::create(target_path)?;

    // Getting data
    for counter in 1..=10000 {
        // let source_path = format!("target/preexplorer/data/sandpiper_poly_{}.txt", counter);
        let source_path = format!("target/preexplorer/data/redneck_poly_{}.txt", counter);
        let input = File::open(source_path)?;
        let buffered = BufReader::new(input);
        for line in buffered.lines().skip(3) {
            write!(output, "{}\n", line?)?;
        }
    }

    // // Checking
    // let input = File::open(target_path)?;
    // let buffered = BufReader::new(input);
    // for line in buffered.lines() {
    //     println!("{}", line?);
    // }

    Ok(())
}
