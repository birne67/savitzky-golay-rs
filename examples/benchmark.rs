 use savitzky_golay::SavitzkyGolayFilter;
use std::time::Instant;

fn main() {
    println!("=== Savitzky-Golay Filter Benchmarks ===\n");

    // Test different data sizes
    let sizes = vec![100, 1_000, 10_000, 100_000];
    
    for size in sizes {
        println!("Benchmarking with {} data points:", size);
        
        // Generate test data
        let data: Vec<f64> = (0..size)
            .map(|i| {
                let x = i as f64 * 0.01;
                x.sin() + 0.5 * (2.0 * x).cos() + 0.1 * (5.0 * x).sin()
            })
            .collect();
        
        // Test different filter configurations
        let configs = vec![
            (5, 2, "5-point quadratic"),
            (7, 3, "7-point cubic"),
            (11, 3, "11-point cubic"),
            (15, 4, "15-point quartic"),
        ];
        
        for (window_size, poly_order, description) in configs {
            let mut filter = SavitzkyGolayFilter::new(window_size, poly_order)
                .expect("Valid parameters");
            
            // Warm up
            let _ = filter.apply(&data[0..std::cmp::min(100, size)]);
            
            // Benchmark smoothing
            let start = Instant::now();
            let _smoothed = filter.apply(&data);
            let smooth_duration = start.elapsed();
            
            // Benchmark first derivative
            let start = Instant::now();
            let _derivative = filter.apply_derivative(&data, 1, 1.0);
            let deriv_duration = start.elapsed();
            
            println!("  {}: smooth={:?}, derivative={:?}", 
                description, smooth_duration, deriv_duration);
        }
        
        println!();
    }
    
    // Memory usage test
    println!("Memory efficiency test:");
    let large_data: Vec<f64> = (0..1_000_000)
        .map(|i| (i as f64 * 0.0001).sin())
        .collect();
    
    println!("  Processing 1M points...");
    let mut filter = SavitzkyGolayFilter::new(11, 3).unwrap();
    
    let start = Instant::now();
    let _result = filter.apply(&large_data);
    let duration = start.elapsed();
    
    println!("  Completed in {:?}", duration);
    println!("  Throughput: {:.0} points/ms", 
        1_000_000.0 / duration.as_millis() as f64);
    
    // Coefficient caching test
    println!("\nCoefficient caching efficiency:");
    let test_data: Vec<f64> = (0..1000).map(|i| i as f64).collect();
    
    // First run (cold cache)
    let mut filter = SavitzkyGolayFilter::new(7, 2).unwrap();
    let start = Instant::now();
    for _ in 0..10 {
        let _ = filter.apply(&test_data);
    }
    let cached_duration = start.elapsed();
    
    println!("  10 runs with caching: {:?}", cached_duration);
    println!("  Average per run: {:?}", cached_duration / 10);
}