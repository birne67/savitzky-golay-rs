use savitzky_golay::*;

fn main() {
    // Debug coefficient calculation
    println!("=== Debugging Savitzky-Golay Implementation ===\n");
    
    // Test linear trend preservation
    println!("Testing linear trend preservation:");
    let linear_data: Vec<f64> = (0..10).map(|i| 2.0 * i as f64 + 3.0).collect();
    println!("Input linear data: {:?}", linear_data);
    
    let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    let filtered = filter.apply(&linear_data);
    println!("Filtered data: {:?}", filtered);
    
    println!("\nDifferences:");
    for (i, (orig, filt)) in linear_data.iter().zip(filtered.iter()).enumerate() {
        println!("  [{}]: {:.6} -> {:.6} (diff: {:.6})", i, orig, filt, orig - filt);
    }
    
    // Test with constant data
    println!("\nTesting constant preservation:");
    let constant_data = vec![5.0; 10];
    println!("Input constant data: {:?}", constant_data);
    
    let filtered_const = filter.apply(&constant_data);
    println!("Filtered data: {:?}", filtered_const);
    
    // Test with unit impulse
    println!("\nTesting unit impulse response:");
    let impulse_data = vec![0.0, 0.0, 1.0, 0.0, 0.0];
    println!("Input impulse data: {:?}", impulse_data);
    
    let mut impulse_filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    let filtered_impulse = impulse_filter.apply(&impulse_data);
    println!("Filtered impulse: {:?}", filtered_impulse);
}