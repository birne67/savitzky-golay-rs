//! # Savitzky-Golay Filter
//!
//! A high-performance implementation of the Savitzky-Golay filter for signal smoothing
//! and numerical differentiation in Rust.
//!
//! The Savitzky-Golay filter is a digital filter that can be applied to a set of 
//! digital data points for the purpose of smoothing the data, that is, to increase 
//! the precision of the data without distorting the signal tendency.
//!
//! ## Features
//!
//! - Fast coefficient calculation using least squares polynomial fitting
//! - Support for different polynomial orders and window sizes
//! - Numerical differentiation capabilities
//! - Boundary handling strategies
//! - Zero-copy operations where possible
//!
//! ## Example
//!
//! ```rust
//! use savitzky_golay::SavitzkyGolayFilter;
//!
//! let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
//! let mut filter = SavitzkyGolayFilter::new(5, 2).expect("Invalid parameters");
//! let smoothed = filter.apply(&data);
//! ```

mod coefficients;
mod filter;
mod error;

pub use filter::{SavitzkyGolayFilter, BoundaryMode, FilterConfig};
pub use error::{SavitzkyGolayError, Result};
pub use coefficients::compute_coefficients;

/// Applies a Savitzky-Golay filter to the input data with default parameters.
///
/// This is a convenience function that creates a filter with window size 5 and 
/// polynomial order 2, which works well for most smoothing applications.
///
/// # Arguments
///
/// * `data` - The input signal data
///
/// # Returns
///
/// A vector containing the filtered data
///
/// # Example
///
/// ```rust
/// use savitzky_golay::smooth;
///
/// let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
/// let smoothed = smooth(&data);
/// ```
pub fn smooth(data: &[f64]) -> Vec<f64> {
    let mut filter = SavitzkyGolayFilter::new(5, 2)
        .expect("Default parameters should be valid");
    filter.apply(data)
}

/// Computes the first derivative of the input data using Savitzky-Golay differentiation.
///
/// # Arguments
///
/// * `data` - The input signal data
/// * `window_size` - The size of the filter window (must be odd)
/// * `poly_order` - The order of the polynomial (must be less than window_size)
///
/// # Returns
///
/// A Result containing the derivative data or an error
pub fn derivative(data: &[f64], window_size: usize, poly_order: usize) -> Result<Vec<f64>> {
    let mut filter = SavitzkyGolayFilter::new(window_size, poly_order)?;
    Ok(filter.apply_derivative(data, 1, 1.0))
}