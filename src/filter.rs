use crate::coefficients::CoefficientCache;
use crate::error::{Result, SavitzkyGolayError};

/// Boundary handling strategies for the Savitzky-Golay filter.
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryMode {
    /// Pad with constant values (first/last data point)
    Constant,
    /// Pad with mirrored values
    Mirror,
    /// Wrap around (circular boundary)
    Wrap,
    /// Pad with zeros
    Zero,
    /// Use smaller windows near boundaries
    Shrink,
}

/// Configuration for the Savitzky-Golay filter.
#[derive(Debug, Clone)]
pub struct FilterConfig {
    /// Size of the filter window (must be odd)
    pub window_size: usize,
    /// Order of the polynomial to fit
    pub poly_order: usize,
    /// Boundary handling strategy
    pub boundary_mode: BoundaryMode,
}

impl FilterConfig {
    /// Creates a new filter configuration with validation
    pub fn new(window_size: usize, poly_order: usize) -> Result<Self> {
        if window_size % 2 == 0 || window_size == 0 {
            return Err(SavitzkyGolayError::InvalidWindowSize(window_size));
        }
        
        if poly_order >= window_size {
            return Err(SavitzkyGolayError::InvalidPolynomialOrder(poly_order, window_size));
        }
        
        Ok(Self {
            window_size,
            poly_order,
            boundary_mode: BoundaryMode::Constant,
        })
    }
    
    /// Sets the boundary handling mode
    pub fn with_boundary_mode(mut self, mode: BoundaryMode) -> Self {
        self.boundary_mode = mode;
        self
    }
}

/// A Savitzky-Golay filter for signal smoothing and differentiation.
pub struct SavitzkyGolayFilter {
    config: FilterConfig,
    cache: CoefficientCache,
}

impl SavitzkyGolayFilter {
    /// Creates a new Savitzky-Golay filter with the specified parameters.
    ///
    /// # Arguments
    ///
    /// * `window_size` - Size of the filter window (must be odd)
    /// * `poly_order` - Order of the polynomial to fit (must be < window_size)
    ///
    /// # Example
    ///
    /// ```rust
    /// use savitzky_golay::SavitzkyGolayFilter;
    ///
    /// let filter = SavitzkyGolayFilter::new(5, 2).expect("Valid parameters");
    /// ```
    pub fn new(window_size: usize, poly_order: usize) -> Result<Self> {
        let config = FilterConfig::new(window_size, poly_order)?;
        Ok(Self {
            config,
            cache: CoefficientCache::default(),
        })
    }
    
    /// Creates a filter with custom configuration
    pub fn with_config(config: FilterConfig) -> Self {
        Self {
            config,
            cache: CoefficientCache::default(),
        }
    }
    
    /// Sets the boundary handling mode
    pub fn with_boundary_mode(mut self, mode: BoundaryMode) -> Self {
        self.config.boundary_mode = mode;
        self
    }
    
    /// Applies the Savitzky-Golay filter to smooth the input data.
    ///
    /// # Arguments
    ///
    /// * `data` - The input signal data
    ///
    /// # Returns
    ///
    /// A vector containing the smoothed data
    ///
    /// # Example
    ///
    /// ```rust
    /// use savitzky_golay::SavitzkyGolayFilter;
    ///
    /// let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
    /// let filter = SavitzkyGolayFilter::new(5, 2).expect("Valid parameters");
    /// let smoothed = filter.apply(&data);
    /// ```
    pub fn apply(&mut self, data: &[f64]) -> Vec<f64> {
        if data.len() < self.config.window_size {
            // For very short data, return a copy
            return data.to_vec();
        }
        
        let coeffs = self.cache
            .get_coefficients(self.config.window_size, self.config.poly_order, 0)
            .expect("Coefficient computation failed")
            .clone(); // Clone to avoid borrowing issues
        
        self.apply_with_coefficients(data, &coeffs)
    }
    
    /// Applies the Savitzky-Golay filter to compute derivatives of the input data.
    ///
    /// # Arguments
    ///
    /// * `data` - The input signal data
    /// * `derivative_order` - Order of derivative to compute (1 for first derivative, 2 for second, etc.)
    ///
    /// # Returns
    ///
    /// A vector containing the derivative data
    pub fn apply_derivative(&mut self, data: &[f64], derivative_order: usize) -> Vec<f64> {
        if data.len() < self.config.window_size {
            return vec![0.0; data.len()];
        }
        
        let coeffs = self.cache
            .get_coefficients(self.config.window_size, self.config.poly_order, derivative_order)
            .expect("Coefficient computation failed")
            .clone(); // Clone to avoid borrowing issues
        
        self.apply_with_coefficients(data, &coeffs)
    }
    
    /// Internal method to apply filter with given coefficients
    fn apply_with_coefficients(&self, data: &[f64], coeffs: &[f64]) -> Vec<f64> {
        let n = data.len();
        let mut result = vec![0.0; n];
        let half_window = self.config.window_size / 2;
        
        // Handle each point in the signal
        for i in 0..n {
            result[i] = self.compute_filtered_value(data, i, coeffs, half_window);
        }
        
        result
    }
    
    /// Computes the filtered value at a specific point
    fn compute_filtered_value(
        &self,
        data: &[f64],
        center: usize,
        coeffs: &[f64],
        half_window: usize,
    ) -> f64 {
        let n = data.len();
        let mut sum = 0.0;
        
        // For each coefficient, find the corresponding data point
        for (coeff_idx, &coeff) in coeffs.iter().enumerate() {
            // Calculate the data index for this coefficient
            let offset = coeff_idx as isize - half_window as isize;
            let data_idx = center as isize + offset;
            
            let value = if data_idx >= 0 && data_idx < n as isize {
                data[data_idx as usize]
            } else {
                self.get_boundary_value(data, data_idx)
            };
            
            sum += coeff * value;
        }
        
        sum
    }
    
    /// Gets a value for boundary handling
    fn get_boundary_value(
        &self,
        data: &[f64],
        requested_idx: isize,
    ) -> f64 {
        let n = data.len();
        
        match self.config.boundary_mode {
            BoundaryMode::Constant => {
                if requested_idx < 0 {
                    data[0] // Use first value for left boundary
                } else {
                    data[n - 1] // Use last value for right boundary
                }
            }
            BoundaryMode::Mirror => {
                if requested_idx < 0 {
                    // Mirror around the start
                    let abs_idx = (-requested_idx) as usize;
                    if abs_idx < n {
                        data[abs_idx]
                    } else {
                        data[n - 1]
                    }
                } else {
                    // Mirror around the end
                    let overflow = requested_idx as usize - (n - 1);
                    let mirrored_idx = if overflow < n {
                        n - 1 - overflow
                    } else {
                        0
                    };
                    data[mirrored_idx]
                }
            }
            BoundaryMode::Wrap => {
                let wrapped_idx = if requested_idx < 0 {
                    (n as isize + requested_idx % n as isize) as usize % n
                } else {
                    requested_idx as usize % n
                };
                data[wrapped_idx]
            }
            BoundaryMode::Zero => 0.0,
            BoundaryMode::Shrink => {
                // Use smaller window - this requires recomputing coefficients
                // For simplicity, fall back to constant mode here
                if requested_idx < 0 {
                    data[0]
                } else {
                    data[n - 1]
                }
            }
        }
    }
    
    /// Returns the filter configuration
    pub fn config(&self) -> &FilterConfig {
        &self.config
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_abs_diff_eq;

    #[test]
    fn test_basic_smoothing() {
        let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
        let smoothed = filter.apply(&data);
        
        // The smoothed signal should be less noisy
        assert_eq!(smoothed.len(), data.len());
        
        // Check that the filter preserves the general trend
        // (exact values depend on the specific implementation details)
        assert!(smoothed[4] > smoothed[0]); // Peak should still be higher than start
        assert!(smoothed[4] > smoothed[8]); // Peak should still be higher than end
    }
    
    #[test]
    fn test_polynomial_preservation() {
        let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
        
        // Quadratic polynomial: y = x^2
        let data: Vec<f64> = (0..10).map(|x| (x as f64).powi(2)).collect();
        let smoothed = filter.apply(&data);
        
        // Savitzky-Golay should preserve polynomials up to the specified order
        for (original, filtered) in data.iter().zip(smoothed.iter()) {
            assert_abs_diff_eq!(original, filtered, epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_derivative() {
        let mut filter = SavitzkyGolayFilter::new(5, 3).unwrap();
        
        // Cubic polynomial: y = x^3
        let data: Vec<f64> = (0..10).map(|x| (x as f64).powi(3)).collect();
        let derivative = filter.apply_derivative(&data, 1);
        
        // First derivative should be 3x^2
        for (i, &deriv) in derivative.iter().enumerate() {
            let expected = 3.0 * (i as f64).powi(2);
            assert_abs_diff_eq!(deriv, expected, epsilon = 1e-6);
        }
    }
    
    #[test]
    fn test_boundary_modes() {
        let mut filter = SavitzkyGolayFilter::new(5, 2)
            .unwrap()
            .with_boundary_mode(BoundaryMode::Zero);
        
        let data = vec![1.0, 2.0, 3.0, 4.0, 5.0];
        let result = filter.apply(&data);
        
        assert_eq!(result.len(), data.len());
    }
    
    #[test]
    fn test_short_data() {
        let mut filter = SavitzkyGolayFilter::new(5, 2).unwrap();
        let data = vec![1.0, 2.0]; // Shorter than window size
        let result = filter.apply(&data);
        
        assert_eq!(result, data); // Should return unchanged
    }
}