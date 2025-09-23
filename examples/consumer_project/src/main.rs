use savitzky_golay::{smooth, SavitzkyGolayFilter, BoundaryMode, read_csv_column};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let data = read_csv_column("data.csv", "A")?;
    println!("{:?}", data);


    // Use convenience function
    let sm = smooth(&data);
    println!("smoothed: {:?}", sm);

    // Use filter directly
    let mut filter = SavitzkyGolayFilter::new(5, 2).expect("valid params");
    filter = filter.with_boundary_mode(BoundaryMode::Mirror);
    let res = filter.apply(&data);
    println!("apply result: {:?}", res);

    // Derivative
    let mut filter_d = SavitzkyGolayFilter::new(5, 3).expect("valid params");
    let deriv = filter_d.apply_derivative(&data, 1, 1.0);
    println!("derivative: {:?}", deriv);

    Ok(())
}
