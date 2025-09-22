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
    // Cache for asymmetric/delta-aware coefficients: key is (hash_of_offsets, poly_order, derivative, delta_bits)
    asymmetric: std::collections::HashMap<(u64, usize, usize, u64), Vec<f64>>,
}

impl CoefficientCache {
    /// Creates a new coefficient cache
    pub fn new() -> Self {
        Self {
            coefficients: std::collections::HashMap::new(),
            asymmetric: std::collections::HashMap::new(),
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

    /// Get asymmetric/delta-aware coefficients, caching them by a hash of offsets, poly order, derivative and delta.
    pub fn get_coefficients_for_offsets_with_delta(
        &mut self,
        offsets: &[isize],
        poly_order: usize,
        derivative: usize,
        delta: f64,
    ) -> Result<&Vec<f64>> {
        // Build a small hash of offsets
        use std::hash::{Hasher, Hash};
        let mut hasher = ahash::AHasher::default();
        offsets.hash(&mut hasher);
        let offsets_hash = hasher.finish();

        let delta_bits = delta.to_bits();
        let key = (offsets_hash, poly_order, derivative, delta_bits);

        if !self.asymmetric.contains_key(&key) {
            let coeffs = compute_coefficients_for_offsets_with_delta(offsets, poly_order, derivative, delta)?;
            self.asymmetric.insert(key, coeffs);
        }

        Ok(&self.asymmetric[&key])
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


/// Applies Savitzky-Golay filter to a signal with specified boundary mode.
///
/// # Arguments
///
/// * `signal` - Input signal as a slice of f64
/// * `coefficients` - Filter coefficients (from `compute_coefficients`)
/// * `mode` - Boundary handling mode
///
/// # Returns
///
/// Filtered signal as a Vec<f64>
/// 
pub fn apply_filter_with_mode(
    signal: &[f64],
    coefficients: &[f64],
    mode: crate::filter::BoundaryMode,
) -> Result<Vec<f64>> {
    if signal.is_empty() || coefficients.is_empty() {
        return Err(SavitzkyGolayError::InvalidInput("Signal or coefficients are empty".to_string()));
    }
    let window_size = coefficients.len();
    let half_window = (window_size - 1) / 2;
    let padded_len = signal.len() + 2 * half_window;
    let mut padded_signal = vec![0.0; padded_len];

    // Pad the signal based on mode
    match mode {
        crate::filter::BoundaryMode::Mirror => {
            // Mirror at boundaries
            for i in 0..half_window {
                padded_signal[i] = signal[half_window - i];
                padded_signal[padded_len - 1 - i] = signal[signal.len() - 1 - (half_window - i)];
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
    crate::filter::BoundaryMode::Constant => {
            // Pad with edge values
            let left_val = signal[0];
            let right_val = signal[signal.len() - 1];
            for i in 0..half_window {
                padded_signal[i] = left_val;
                padded_signal[padded_len - 1 - i] = right_val;
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
    crate::filter::BoundaryMode::Nearest => {
            // Repeat nearest values (similar to constant but per position)
            for i in 0..half_window {
                padded_signal[i] = signal[0];
                padded_signal[padded_len - 1 - i] = signal[signal.len() - 1];
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
    crate::filter::BoundaryMode::Wrap => {
            // Wrap around
            for i in 0..half_window {
                padded_signal[i] = signal[signal.len() - half_window + i];
                padded_signal[padded_len - 1 - i] = signal[i];
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
    crate::filter::BoundaryMode::Interp => {
            // Linear interpolation/extrapolation
            if signal.len() >= 2 {
                let left_slope = signal[1] - signal[0];
                let right_slope = signal[signal.len() - 1] - signal[signal.len() - 2];
                for i in 0..half_window {
                    padded_signal[i] = signal[0] - left_slope * (half_window - i) as f64;
                    padded_signal[padded_len - 1 - i] = signal[signal.len() - 1] + right_slope * (i + 1) as f64;
                }
            } else {
                // Fallback to constant if not enough points
                let val = *signal.get(0).unwrap_or(&0.0);
                padded_signal.fill(val);
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
        crate::filter::BoundaryMode::Zero => {
            // Pad with zeros
            for i in 0..padded_len {
                padded_signal[i] = 0.0;
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
        crate::filter::BoundaryMode::Shrink => {
            // Fall back to constant padding for now
            let left_val = signal.get(0).copied().unwrap_or(0.0);
            let right_val = signal.last().copied().unwrap_or(0.0);
            for i in 0..half_window {
                padded_signal[i] = left_val;
                padded_signal[padded_len - 1 - i] = right_val;
            }
            padded_signal[half_window..half_window + signal.len()].copy_from_slice(signal);
        }
    }

    // Convolve with coefficients
    let mut filtered = vec![0.0; signal.len()];
    for i in 0..signal.len() {
        for j in 0..window_size {
            filtered[i] += padded_signal[i + j] * coefficients[j];
        }
    }

    Ok(filtered)
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

    #[test]
    fn test_apply_filter_mirror_mode() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0];
    let filtered = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Mirror).unwrap();
        
        // Check that filtered signal has same length as input
        assert_eq!(filtered.len(), signal.len());
        // Basic sanity check - filtered values should be reasonable
        for &val in &filtered {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_apply_filter_constant_mode() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let filtered = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Constant).unwrap();
        
        assert_eq!(filtered.len(), signal.len());
        for &val in &filtered {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_apply_filter_nearest_mode() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let filtered = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Nearest).unwrap();
        
        assert_eq!(filtered.len(), signal.len());
        for &val in &filtered {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_apply_filter_wrap_mode() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let filtered = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Wrap).unwrap();
        
        assert_eq!(filtered.len(), signal.len());
        for &val in &filtered {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_apply_filter_interp_mode() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let filtered = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Interp).unwrap();
        
        assert_eq!(filtered.len(), signal.len());
        for &val in &filtered {
            assert!(val.is_finite());
        }
    }

    #[test]
    fn test_apply_filter_empty_signal() {
        let coeffs = compute_coefficients(5, 2, 0).unwrap();
        let signal = vec![];
    let result = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Mirror);
        assert!(result.is_err());
    }

    #[test]
    fn test_apply_filter_empty_coeffs() {
        let coeffs = vec![];
        let signal = vec![1.0, 2.0, 3.0];
    let result = apply_filter_with_mode(&signal, &coeffs, crate::filter::BoundaryMode::Mirror);
        assert!(result.is_err());
    }
}

// Non-delta offsets variant removed; use `compute_coefficients_for_offsets_with_delta` for
// asymmetric coefficient computation (pass delta=1.0 for unit spacing).

/// Compute coefficients taking the sample spacing `delta` into account.
/// Using `delta` here builds the Vandermonde in the real x-units so the returned
/// filter coefficients directly approximate the derivative of order `derivative`
/// in physical units (no later scaling required).
pub fn compute_coefficients_with_delta(
    window_size: usize,
    poly_order: usize,
    derivative: usize,
    delta: f64,
) -> Result<Vec<f64>> {
    // Validate parameters (reuse existing checks)
    if window_size % 2 == 0 || window_size == 0 {
        return Err(SavitzkyGolayError::InvalidWindowSize(window_size));
    }
    if poly_order >= window_size {
        return Err(SavitzkyGolayError::InvalidPolynomialOrder(poly_order, window_size));
    }

    let half_window = (window_size - 1) / 2;
    let num_points = window_size;

    // Build Vandermonde in x-units (x = offset * delta)
    let mut vandermonde = DMatrix::<f64>::zeros(num_points, poly_order + 1);
    for i in 0..num_points {
        let x = ((i as isize) - (half_window as isize)) as f64 * delta;
        for j in 0..=poly_order {
            vandermonde[(i, j)] = x.powi(j as i32);
        }
    }

    let ata = vandermonde.transpose() * &vandermonde;
    let mut rhs = DVector::<f64>::zeros(poly_order + 1);

    if derivative <= poly_order {
        let factorial = (1..=derivative).fold(1.0, |acc, x| acc * x as f64);
        rhs[derivative] = factorial;
    } else {
        return Ok(vec![0.0; window_size]);
    }

    let coeffs_poly = ata.lu().solve(&rhs)
        .ok_or_else(|| SavitzkyGolayError::ComputationError(
            "Failed to solve least squares system".to_string()
        ))?;

    let mut filter_coeffs = vec![0.0; window_size];
    for i in 0..window_size {
        let x = ((i as isize) - (half_window as isize)) as f64 * delta;
        for j in 0..=poly_order {
            filter_coeffs[i] += coeffs_poly[j] * x.powi(j as i32);
        }
    }

    Ok(filter_coeffs)
}

/// As above but accepts `delta` to express offsets in physical units.
pub fn compute_coefficients_for_offsets_with_delta(
    offsets: &[isize],
    poly_order: usize,
    derivative: usize,
    delta: f64,
) -> Result<Vec<f64>> {
    let window_size = offsets.len();
    if window_size == 0 {
        return Err(SavitzkyGolayError::InvalidWindowSize(0));
    }

    if poly_order >= window_size {
        return Err(SavitzkyGolayError::InvalidPolynomialOrder(poly_order, window_size));
    }

    let mut vandermonde = DMatrix::<f64>::zeros(window_size, poly_order + 1);
    for (i, &off) in offsets.iter().enumerate() {
        let x = off as f64 * delta;
        for j in 0..=poly_order {
            vandermonde[(i, j)] = x.powi(j as i32);
        }
    }

    let ata = vandermonde.transpose() * &vandermonde;
    let mut rhs = DVector::<f64>::zeros(poly_order + 1);

    if derivative <= poly_order {
        let factorial = (1..=derivative).fold(1.0, |acc, x| acc * x as f64);
        rhs[derivative] = factorial;
    } else {
        return Ok(vec![0.0; window_size]);
    }

    let coeffs_poly = ata.lu().solve(&rhs)
        .ok_or_else(|| SavitzkyGolayError::ComputationError(
            "Failed to solve least squares system (offsets)".to_string()
        ))?;

    let mut filter_coeffs = vec![0.0; window_size];
    for i in 0..window_size {
        let x = offsets[i] as f64 * delta;
        for j in 0..=poly_order {
            filter_coeffs[i] += coeffs_poly[j] * x.powi(j as i32);
        }
    }

    Ok(filter_coeffs)
}

/// Backwards-compatible wrapper: compute coefficients for arbitrary offsets using unit spacing.
/// This preserves the old API while keeping the canonical implementation in
/// `compute_coefficients_for_offsets_with_delta`.
pub fn compute_coefficients_for_offsets(
    offsets: &[isize],
    poly_order: usize,
    derivative: usize,
) -> Result<Vec<f64>> {
    compute_coefficients_for_offsets_with_delta(offsets, poly_order, derivative, 1.0)
}