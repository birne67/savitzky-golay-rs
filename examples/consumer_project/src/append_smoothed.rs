
use std::error::Error;
use std::fs::File;
use std::path::Path;
use savitzky_golay::{SavitzkyGolayFilter, BoundaryMode, smooth, derivative, read_csv_column};

fn main() -> Result<(), Box<dyn Error>> {
    let input_path = "data.csv";
    let output_path = "data.csv"; // Overwrite the same file
    let column_name = "A"; // Change as needed
    let smoothed_col_name = "smoothed";

    // Read the whole CSV into memory
    let mut rdr = csv::Reader::from_path(input_path)?;
    let headers = rdr.headers()?.clone();
    let mut records: Vec<csv::StringRecord> = rdr.records().map(|r| r.unwrap()).collect();

    // Read the target column as f64
    let data = read_csv_column(input_path, column_name)?;

    let mut nearest_filter = SavitzkyGolayFilter::new(41, 5)?
        .with_boundary_mode(BoundaryMode::Nearest);
    let smoothed = nearest_filter.apply(&data);


    // Add the new column to each record
    for (i, record) in records.iter_mut().enumerate() {
        // If smoothed value exists, append it; else, append empty
        let value = smoothed.get(i).map(|v: &f64| v.to_string()).unwrap_or_default();
        record.push_field(&value);
    }

    // Write out the new CSV
    let mut wtr = csv::Writer::from_path(output_path)?;
    let mut new_headers = headers.clone();
    new_headers.push_field(smoothed_col_name);
    wtr.write_record(&new_headers)?;
    for record in records {
        wtr.write_record(&record)?;
    }
    wtr.flush()?;
    println!("Smoothed column appended to {}!", output_path);
    Ok(())
}
