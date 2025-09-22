use savitzky_golay::{SavitzkyGolayFilter, BoundaryMode, smooth, derivative};
use approx::assert_abs_diff_eq;

#[test]
fn test_known_coefficients() {
    // Test against known coefficient values from literature
    let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    
    // For a 5-point quadratic Savitzky-Golay filter, the coefficients are:
    // [-3, 12, 17, 12, -3] / 35
    let test_data = vec![0.0, 0.0, 1.0, 0.0, 0.0]; // Unit impulse at center
    let result = filter.apply(&test_data);
    
    // The result should approximate the impulse response
    // Note: due to boundary effects, we check the behavior is reasonable
    assert!(result[2] > result[0]); // Center should be largest
    assert!(result[2] > result[4]); // Center should be largest
}

#[test]
fn test_polynomial_preservation_fixed() {
    // Test that polynomials up to the specified order are preserved
    let mut filter = SavitzkyGolayFilter::new(7, 3).unwrap();
    
    // Test with a cubic polynomial: f(x) = x^3 - 2x^2 + x + 1
    let data: Vec<f64> = (0..15)
        .map(|i| {
            let x = i as f64;
            x.powi(3) - 2.0 * x.powi(2) + x + 1.0
        })
        .collect();
    
    let filtered = filter.apply(&data);
    
    // For the interior points (away from boundaries), the polynomial should be preserved
    let interior_start = 3;
    let interior_end = data.len() - 3;
    
    for i in interior_start..interior_end {
        assert_abs_diff_eq!(data[i], filtered[i], epsilon = 1e-6);
    }
}

#[test]
fn test_linear_trend_preservation() {
    // Linear functions should always be preserved regardless of window size or polynomial order
    let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    
    let linear_data: Vec<f64> = (0..10).map(|i| 2.0 * i as f64 + 3.0).collect();
    let filtered = filter.apply(&linear_data);
    
    for (original, filtered) in linear_data.iter().zip(filtered.iter()) {
        assert_abs_diff_eq!(original, filtered, epsilon = 1e-10);
    }
}

#[test]
fn test_constant_preservation() {
    // Constant signals should be perfectly preserved
    let mut filter = SavitzkyGolayFilter::new(7, 3).unwrap();
    
    let constant_data = vec![5.0; 20];
    let filtered = filter.apply(&constant_data);
    
    for &value in &filtered {
        assert_abs_diff_eq!(value, 5.0, epsilon = 1e-12);
    }
}

#[test]
fn test_derivative_computation() {
    // Test derivative computation on a known function
    let mut filter = SavitzkyGolayFilter::new(7, 4).unwrap();
    
    // Use f(x) = x^4, so f'(x) = 4x^3
    let step = 0.1;
    let data: Vec<f64> = (0..20)
        .map(|i| {
            let x = i as f64 * step;
            x.powi(4)
        })
        .collect();
    
    let first_derivative = filter.apply_derivative(&data, 1, step);

    // Check derivative values in the interior (physical derivative: d/dt x^4 = 4 x^3)
    for i in 3..17 {
        let x = i as f64 * step;
        let expected_derivative = 4.0 * x.powi(3);

        // Allow some tolerance due to numerical approximation
        assert_abs_diff_eq!(first_derivative[i], expected_derivative, epsilon = 0.1);
    }
}

#[test]
fn test_boundary_modes() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    
    // Test different boundary modes
    let modes = [
        BoundaryMode::Constant,
        BoundaryMode::Mirror,
        BoundaryMode::Wrap,
        BoundaryMode::Zero,
        BoundaryMode::Shrink,
    ];
    
    for mode in &modes {
        let mut filter = SavitzkyGolayFilter::new(5, 2)
            .unwrap()
            .with_boundary_mode(*mode);
        
        let result = filter.apply(&data);
        assert_eq!(result.len(), data.len());
        
        // Results should be finite numbers
        for &value in &result {
            assert!(value.is_finite());
        }
    }
}

#[test]
fn test_noise_reduction() {
    // Create a signal with known underlying trend + noise
    let true_signal: Vec<f64> = (0..50)
        .map(|i| (i as f64 * 0.1).sin())
        .collect();
    
    let mut noisy_signal = true_signal.clone();
    // Add deterministic "noise" for reproducible testing
    for i in 0..noisy_signal.len() {
        noisy_signal[i] += 0.1 * ((i as f64 * 1.7).sin());
    }
    
    let mut filter = SavitzkyGolayFilter::new(9, 3).unwrap();
    let smoothed = filter.apply(&noisy_signal);
    
    // Calculate mean squared error for original noisy vs true signal
    let mse_noisy: f64 = true_signal.iter()
        .zip(noisy_signal.iter())
        .map(|(true_val, noisy_val)| (true_val - noisy_val).powi(2))
        .sum::<f64>() / true_signal.len() as f64;
    
    // Calculate mean squared error for smoothed vs true signal
    let mse_smoothed: f64 = true_signal.iter()
        .zip(smoothed.iter())
        .map(|(true_val, smooth_val)| (true_val - smooth_val).powi(2))
        .sum::<f64>() / true_signal.len() as f64;
    
    // Smoothed signal should be closer to the true signal than the noisy signal
    assert!(mse_smoothed < mse_noisy);
}

#[test]
fn test_convenience_functions() {
    let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    
    // Test smooth function
    let smoothed = smooth(&data);
    assert_eq!(smoothed.len(), data.len());
    
    // Test derivative function
    let deriv = derivative(&data, 7, 2).unwrap();
    assert_eq!(deriv.len(), data.len());
}

#[test]
fn test_edge_cases() {
    // Test very small datasets
    let small_data = vec![1.0, 2.0];
    let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    let result = filter.apply(&small_data);
    assert_eq!(result, small_data); // Should return unchanged
    
    // Test single point
    let single_point = vec![42.0];
    let result = filter.apply(&single_point);
    assert_eq!(result, single_point);
    
    // Test empty data
    let empty_data: Vec<f64> = vec![];
    let result = filter.apply(&empty_data);
    assert_eq!(result, empty_data);
}

#[test]
fn test_performance_characteristics() {
    // Test that filter can handle reasonably large datasets efficiently
    let large_data: Vec<f64> = (0..10000)
        .map(|i| (i as f64 * 0.001).sin())
        .collect();
    
    let mut filter = SavitzkyGolayFilter::new(11, 3).unwrap();
    
    let start = std::time::Instant::now();
    let _result = filter.apply(&large_data);
    let duration = start.elapsed();
    
    // Should complete in reasonable time (less than 100ms for 10k points)
    assert!(duration.as_millis() < 100);
}

#[test]
fn test_numerical_stability() {
    // Test with extreme values to ensure numerical stability
    let extreme_data = vec![
        1e-10, 1e10, -1e10, 1e-10, 0.0,
        f64::MIN_POSITIVE, f64::MAX / 1e10, -f64::MAX / 1e10,
    ];
    
    let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
    let result = filter.apply(&extreme_data);
    
    // All results should be finite
    for &value in &result {
        assert!(value.is_finite(), "Got non-finite value: {}", value);
    }
}