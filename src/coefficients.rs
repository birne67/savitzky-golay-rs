use nalgebra::{DMatrix, DVector};
use crate::error::{Result, SavitzkyGolayError};

/// Computes Savitzky-Golay filter coefficients using least squares polynomial fitting.
///
/// The Savitzky-Golay filter works by fitting a polynomial of degree `poly_order` 
/// to a moving window of `window_size` data points, then using the polynomial to 
/// estimate the value (or derivative) at the center point.
///
/// # Arguments
///
/// * `window_size` - Size of the moving window (must be odd)
/// * `poly_order` - Degree of the polynomial to fit (must be < window_size)
/// * `derivative` - Order of derivative to compute (0 for smoothing, 1 for first derivative, etc.)
///
/// # Returns
///
/// A vector of filter coefficients that can be convolved with the signal
pub fn compute_coefficients(
    window_size: usize,
    poly_order: usize,
    derivative: usize,
) -> Result<Vec<f64>> {
    // Validate parameters
    if window_size % 2 == 0 || window_size == 0 {
        return Err(SavitzkyGolayError::InvalidWindowSize(window_size));
    }
    
    if poly_order >= window_size {
        return Err(SavitzkyGolayError::InvalidPolynomialOrder(poly_order, window_size));
    }

    let half_window = (window_size - 1) / 2;
    let num_points = window_size;
    
    // Create the Vandermonde matrix for polynomial fitting
    // Each row corresponds to a point in the window, each column to a polynomial power
    let mut vandermonde = DMatrix::<f64>::zeros(num_points, poly_order + 1);
    
    for i in 0..num_points {
        let x = (i as f64) - (half_window as f64); // Center the window around 0
        for j in 0..=poly_order {
            vandermonde[(i, j)] = x.powi(j as i32);
        }
    }
    
    // Solve the least squares problem: A^T A c = A^T e_k
    // where e_k is the k-th unit vector (k = derivative order)
    let ata = vandermonde.transpose() * &vandermonde;
    let mut rhs = DVector::<f64>::zeros(poly_order + 1);
    
    // Set up right-hand side for the desired derivative
    if derivative <= poly_order {
        // The coefficient for x^derivative term, scaled by derivative factorial
        let factorial = (1..=derivative).fold(1.0, |acc, x| acc * x as f64);
        rhs[derivative] = factorial;
    } else {
        // Derivative order higher than polynomial order results in zero
        return Ok(vec![0.0; window_size]);
    }
    
    // Solve the system
    let coeffs_poly = ata.lu().solve(&rhs)
        .ok_or_else(|| SavitzkyGolayError::ComputationError(
            "Failed to solve least squares system".to_string()
        ))?;
    
    // Convert polynomial coefficients to filter coefficients
    let mut filter_coeffs = vec![0.0; window_size];
    
    for i in 0..window_size {
        let x = (i as f64) - (half_window as f64);
        for j in 0..=poly_order {
            filter_coeffs[i] += coeffs_poly[j] * x.powi(j as i32);
        }
    }
    
    Ok(filter_coeffs)
}

/// Precomputed coefficients for common Savitzky-Golay filter configurations.
/// This provides faster access to frequently used coefficient sets.
pub struct CoefficientCache {
    coefficients: std::collections::HashMap<(usize, usize, usize), Vec<f64>>,
}

impl CoefficientCache {
    /// Creates a new coefficient cache
    pub fn new() -> Self {
        Self {
            coefficients: std::collections::HashMap::new(),
        }
    }
    
    /// Gets coefficients from cache or computes them if not cached
    pub fn get_coefficients(
        &mut self,
        window_size: usize,
        poly_order: usize,
        derivative: usize,
    ) -> Result<&Vec<f64>> {
        let key = (window_size, poly_order, derivative);
        
        if !self.coefficients.contains_key(&key) {
            let coeffs = compute_coefficients(window_size, poly_order, derivative)?;
            self.coefficients.insert(key, coeffs);
        }
        
        Ok(&self.coefficients[&key])
    }
    
    /// Precomputes coefficients for common configurations
    pub fn precompute_common(&mut self) -> Result<()> {
        let common_configs = [
            (5, 2, 0),   // 5-point quadratic smoothing
            (7, 2, 0),   // 7-point quadratic smoothing
            (9, 2, 0),   // 9-point quadratic smoothing
            (5, 3, 0),   // 5-point cubic smoothing
            (7, 3, 0),   // 7-point cubic smoothing
            (5, 2, 1),   // 5-point quadratic first derivative
            (7, 2, 1),   // 7-point quadratic first derivative
            (5, 2, 2),   // 5-point quadratic second derivative
        ];
        
        for (window_size, poly_order, derivative) in &common_configs {
            self.get_coefficients(*window_size, *poly_order, *derivative)?;
        }
        
        Ok(())
    }
}

impl Default for CoefficientCache {
    fn default() -> Self {
        let mut cache = Self::new();
        let _ = cache.precompute_common(); // Ignore errors during default initialization
        cache
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_5_point_quadratic_smoothing() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        // Known coefficients for 5-point quadratic smoothing
        let expected = vec![-3.0/35.0, 12.0/35.0, 17.0/35.0, 12.0/35.0, -3.0/35.0];
        
        for (actual, expected) in coeffs.iter().zip(expected.iter()) {
            assert_abs_diff_eq!(actual, expected, epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_invalid_parameters() {
        assert!(compute_coefficients(4, 2, 0).is_err()); // Even window size
        assert!(compute_coefficients(5, 5, 0).is_err()); // Poly order >= window size
        assert!(compute_coefficients(0, 2, 0).is_err()); // Zero window size
    }
    
    #[test]
    fn test_coefficient_cache() {
        let mut cache = CoefficientCache::new();
        
        let _coeffs1 = cache.get_coefficients(5, 2, 0).unwrap();
        let _coeffs2 = cache.get_coefficients(5, 2, 0).unwrap();
        
        // Verify that the cache contains the entry
        assert!(cache.coefficients.contains_key(&(5, 2, 0)));
    }
}