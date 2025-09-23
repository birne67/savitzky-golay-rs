use std::error::Error;
use std::fs::File;
use std::path::Path;

/// Reads a column of f64 values from a CSV file by column name, skipping invalid/missing values.
/// Returns a Vec<f64> with the parsed values.
pub fn read_csv_column<P: AsRef<Path>>(
    path: P,
    column: &str,
) -> Result<Vec<f64>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::Reader::from_reader(file);
    let mut values = Vec::new();

    // Read headers once and get the column index
    let headers = rdr.headers()?.clone();
    let col_index = headers.iter().position(|h| h == column).unwrap_or(0);

    for result in rdr.records() {
        let record = result?;
        if let Some(field) = record.get(col_index) {
            if let Ok(val) = field.parse::<f64>() {
                values.push(val);
            }
        }
    }
    Ok(values)
}

/// Reads a column of f64 values from a CSV file by column index, skipping invalid/missing values.
pub fn read_csv_column_by_index<P: AsRef<Path>>(
    path: P,
    col_index: usize,
) -> Result<Vec<f64>, Box<dyn Error>> {
    let file = File::open(path)?;
    let mut rdr = csv::Reader::from_reader(file);
    let mut values = Vec::new();

    for result in rdr.records() {
        let record = result?;
        if let Some(field) = record.get(col_index) {
            if let Ok(val) = field.parse::<f64>() {
                values.push(val);
            }
        }
    }
    Ok(values)
}
