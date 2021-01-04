use std::fs::File;
use std::io::Error;

fn main() -> Result<(), Error> {
    // if true {
    //     let target_path = "target/preexplorer/data/redneck_poly_all.txt";
    //     let mut output = File::create(target_path)?;

    //     // Getting data
    //     for counter in 1..=10 {
    //         let source_path = format!("target/preexplorer/data/redneck_poly_{}.txt", counter);
    //         let input = File::open(source_path)?;
    //         let buffered = BufReader::new(input);
    //         for line in buffered.lines().skip(3) {
    //             write!(output, "{}\n", line?)?;
    //         }
    //     }
    // }

    if true {
        if true {
            let target_path = "all_redneck_poly.csv";
            let output_file = File::create(target_path)?;
            let mut writer = csv::Writer::from_writer(output_file);

            // Getting data
            for counter in 1..=16_200 {
                collect_record(&mut writer, "redneck", counter)?
            }
            writer.flush()?;
        }
        if true {
            let target_path = "all_sandpiper_poly.csv";
            let output_file = File::create(target_path)?;
            let mut writer = csv::Writer::from_writer(output_file);

            // Getting data
            for counter in 1..=16_200 {
                collect_record(&mut writer, "sandpiper", counter)?
            }
            writer.flush()?;
        }
    
    }

    Ok(())
}

fn collect_record(writer: &mut csv::Writer<File>, bird: &str, counter: usize) -> Result<(), Error> {
    let source_path = format!("{}_poly_{}.csv", bird, counter);
    let mut rdr = csv::ReaderBuilder::new()
        .has_headers(false)
        .from_path(source_path)?;

    if let Some(result) = rdr.records().next() {
        let record = result?;
        writer.write_record(&record)?;
    }
    Ok(())
}