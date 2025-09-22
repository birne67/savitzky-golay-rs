//! Example usage of the Savitzky-Golay filter crate

use savitzky_golay::{SavitzkyGolayFilter, BoundaryMode, smooth, derivative};

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== Savitzky-Golay Filter Examples ===\n");

    // Create some noisy test data
    let clean_signal: Vec<f64> = (0..20)
        .map(|i| {
            let x = i as f64 * 0.1;
            (2.0 * std::f64::consts::PI * x).sin() + 0.5 * (4.0 * std::f64::consts::PI * x).cos()
        })
        .collect();

    // Add some noise
    let mut noisy_signal = clean_signal.clone();
    for i in 0..noisy_signal.len() {
        if i % 3 == 0 {
            noisy_signal[i] += 0.3 * (i as f64 % 2.0 - 0.5);
        }
    }

    println!("Original noisy signal:");
    print_signal(&noisy_signal);

    // Example 1: Basic smoothing with convenience function
    println!("\n1. Basic smoothing (window=5, poly_order=2):");
    let smoothed_basic = smooth(&noisy_signal);
    print_signal(&smoothed_basic);

    // Example 2: Custom filter configuration
    println!("\n2. Custom smoothing (window=7, poly_order=3):");
    let mut custom_filter = SavitzkyGolayFilter::new(7, 3)?;
    let smoothed_custom = custom_filter.apply(&noisy_signal);
    print_signal(&smoothed_custom);

    // Example 3: Different boundary modes
    println!("\n3. Mirror boundary mode:");
    let mut mirror_filter = SavitzkyGolayFilter::new(5, 2)?
        .with_boundary_mode(BoundaryMode::Mirror);
    let smoothed_mirror = mirror_filter.apply(&noisy_signal);
    print_signal(&smoothed_mirror);

    // Example 3b: Nearest boundary mode
    println!("\n3b. Nearest boundary mode:");
    let mut nearest_filter = SavitzkyGolayFilter::new(5, 2)?
        .with_boundary_mode(BoundaryMode::Nearest);
    let smoothed_nearest = nearest_filter.apply(&noisy_signal);
    print_signal(&smoothed_nearest);

    // Example 4: First derivative
    println!("\n4. First derivative:");
    let first_deriv = derivative(&clean_signal, 7, 2)?;
    print_signal(&first_deriv);

    // Example 5: Second derivative
    println!("\n5. Second derivative:");
    let mut deriv_filter = SavitzkyGolayFilter::new(7, 3)?;
    let second_deriv = deriv_filter.apply_derivative(&clean_signal, 2, 1.0);
    print_signal(&second_deriv);

    // Example 6: Performance comparison
    println!("\n6. Performance test with large dataset:");
    let large_data: Vec<f64> = (0..10000)
        .map(|i| (i as f64 * 0.001).sin() + 0.1 * (i as f64 * 0.01).cos())
        .collect();

    let start = std::time::Instant::now();
    let mut perf_filter = SavitzkyGolayFilter::new(11, 3)?;
    let _smoothed_large = perf_filter.apply(&large_data);
    let duration = start.elapsed();
    println!("Processed {} points in {:?}", large_data.len(), duration);

    Ok(())
}

fn print_signal(signal: &[f64]) {
    for (i, &value) in signal.iter().enumerate() {
        print!("{:6.3}", value);
        if i > 0 && (i + 1) % 8 == 0 {
            println!();
        }
    }
    if signal.len() % 8 != 0 {
        println!();
    }
}