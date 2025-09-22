use std::fmt;

/// Error types for Savitzky-Golay filter operations
#[derive(Debug, Clone, PartialEq)]
pub enum SavitzkyGolayError {
    /// Window size must be odd and greater than polynomial order
    InvalidWindowSize(usize),
    /// Polynomial order must be less than window size
    InvalidPolynomialOrder(usize, usize),
    /// Input data is too short for the specified window size
    InsufficientData(usize, usize),
    /// Mathematical computation error (e.g., singular matrix)
    ComputationError(String),
}

impl fmt::Display for SavitzkyGolayError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            SavitzkyGolayError::InvalidWindowSize(size) => {
                write!(f, "Invalid window size: {}. Window size must be odd and at least 1", size)
            }
            SavitzkyGolayError::InvalidPolynomialOrder(poly_order, window_size) => {
                write!(
                    f,
                    "Invalid polynomial order: {}. Must be less than window size ({})",
                    poly_order, window_size
                )
            }
            SavitzkyGolayError::InsufficientData(data_len, window_size) => {
                write!(
                    f,
                    "Insufficient data: {} points. Need at least {} points for window size {}",
                    data_len, window_size, window_size
                )
            }
            SavitzkyGolayError::ComputationError(msg) => {
                write!(f, "Computation error: {}", msg)
            }
        }
    }
}

impl std::error::Error for SavitzkyGolayError {}

/// Result type for Savitzky-Golay operations
pub type Result<T> = std::result::Result<T, SavitzkyGolayError>;