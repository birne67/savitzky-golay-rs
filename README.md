# Savitzky-Golay Filter for Rust

[![Crates.io](https://img.shields.io/crates/v/savitzky_golay.svg)](https://crates.io/crates/savitzky_golay)
[![Documentation](https://docs.rs/savitzky_golay/badge.svg)](https://docs.rs/savitzky_golay)
[![License: MIT OR Apache-2.0](https://img.shields.io/crates/l/savitzky_golay.svg)](#license)

A high-performance implementation of the Savitzky-Golay filter for signal smoothing and numerical differentiation in Rust.

## What is a Savitzky-Golay Filter?

The Savitzky-Golay filter is a digital filter that can be applied to a set of digital data points for the purpose of smoothing the data, increasing the precision of the data without distorting the signal tendency. It achieves this by fitting successive sub-sets of adjacent data points with a low-degree polynomial by the method of linear least squares.

## Features

- 🚀 **High Performance**: Fast coefficient calculation and efficient filtering
- 🔧 **Flexible Configuration**: Customizable window sizes and polynomial orders
- 📊 **Derivative Support**: Compute numerical derivatives of any order
- 🛡️ **Boundary Handling**: Multiple strategies for handling signal boundaries
- ✅ **Well Tested**: Comprehensive test suite with known mathematical results
- 📚 **Zero Dependencies**: Uses only `nalgebra` for linear algebra operations

## Quick Start

Add this to your `Cargo.toml`:

```toml
[dependencies]
savitzky_golay = "0.1"
```

### Basic Usage

```rust
use savitzky_golay::{SavitzkyGolayFilter, smooth, derivative};

// Quick smoothing with default parameters (window=5, poly_order=2)
let data = vec![1.0, 2.0, 3.0, 4.0, 5.0, 4.0, 3.0, 2.0, 1.0];
let smoothed = smooth(&data);

// Custom filter configuration
let mut filter = SavitzkyGolayFilter::new(7, 3)?;
let smoothed = filter.apply(&data);

// Compute first derivative
let first_deriv = derivative(&data, 7, 2)?;

// Compute second derivative
let mut filter = SavitzkyGolayFilter::new(7, 3)?;
let second_deriv = filter.apply_derivative(&data, 2, 1.0);
```

### Advanced Configuration

```rust
use savitzky_golay::{SavitzkyGolayFilter, BoundaryMode};

// Configure boundary handling
let mut filter = SavitzkyGolayFilter::new(5, 2)?
    .with_boundary_mode(BoundaryMode::Mirror);

let smoothed = filter.apply(&data);
```

## Boundary Handling

The filter supports several boundary handling strategies:

- **`Constant`**: Pad with the first/last values (default)
- **`Mirror`**: Reflect the signal at boundaries
- **`Wrap`**: Treat the signal as periodic
- **`Zero`**: Pad with zeros
- **`Shrink`**: Use smaller windows near boundaries

## Performance

The implementation is optimized for performance:

```rust
// Example: Processing 10,000 points typically takes ~3ms
let large_data: Vec<f64> = (0..10000)
    .map(|i| (i as f64 * 0.001).sin())
    .collect();

let mut filter = SavitzkyGolayFilter::new(11, 3)?;
let smoothed = filter.apply(&large_data); // ~3ms
```

## Examples

### Noise Reduction

```rust
use savitzky_golay::SavitzkyGolayFilter;

// Original noisy signal
let noisy_signal = vec![1.1, 1.9, 3.2, 3.8, 5.1, 4.8, 2.9, 2.1, 0.9];

// Apply smoothing
let mut filter = SavitzkyGolayFilter::new(5, 2)?;
let smoothed = filter.apply(&noisy_signal);

println!("Original: {:?}", noisy_signal);
println!("Smoothed: {:?}", smoothed);
```

### Signal Differentiation

```rust
use savitzky_golay::SavitzkyGolayFilter;

// Generate a quadratic signal: y = x²
let signal: Vec<f64> = (0..20)
    .map(|i| (i as f64 * 0.1).powi(2))
    .collect();

// Compute first derivative (should be ≈ 2x)
let mut filter = SavitzkyGolayFilter::new(7, 3)?;
let first_deriv = filter.apply_derivative(&signal, 1, 1.0);

// Compute second derivative (should be ≈ 2)
let second_deriv = filter.apply_derivative(&signal, 2, 1.0);
```

### Multiple Filter Configurations

```rust
use savitzky_golay::{SavitzkyGolayFilter, BoundaryMode};

let data = vec![/* your data */];

// Light smoothing
let mut light_filter = SavitzkyGolayFilter::new(5, 2)?;
let light_smoothed = light_filter.apply(&data);

// Heavy smoothing  
let mut heavy_filter = SavitzkyGolayFilter::new(11, 3)?;
let heavy_smoothed = heavy_filter.apply(&data);

// Custom boundary handling
let mut mirror_filter = SavitzkyGolayFilter::new(7, 2)?
    .with_boundary_mode(BoundaryMode::Mirror);
let mirror_smoothed = mirror_filter.apply(&data);
```

## Mathematical Background

The Savitzky-Golay filter works by:

1. **Local Polynomial Fitting**: For each point, fit a polynomial of degree `p` to `2m+1` surrounding points
2. **Least Squares Solution**: Solve the overdetermined system using least squares
3. **Point Estimation**: Use the fitted polynomial to estimate the value (or derivative) at the center point

The filter coefficients are computed by solving:

```
A^T A c = A^T e_k
```

Where:
- `A` is the Vandermonde matrix of the local window
- `c` are the filter coefficients
- `e_k` is the k-th unit vector (k = derivative order)

## Algorithm Properties

- **Polynomial Preservation**: Preserves polynomials up to the specified order exactly
- **Optimal Smoothing**: Minimizes least squares error for the given polynomial order
- **Local Operation**: Each output point depends only on nearby input points
- **Linear Phase**: Symmetric filters have linear phase response

## When to Use Savitzky-Golay Filters

✅ **Good for:**
- Smoothing signals while preserving features
- Computing derivatives of noisy signals
- Data where you know the underlying trend is polynomial
- Real-time applications (due to computational efficiency)

⚠️ **Consider alternatives for:**
- Signals with sharp discontinuities (may cause ringing)
- Very noisy data (may need pre-filtering)
- Non-polynomial underlying functions (may introduce bias)

## Error Handling

The library provides comprehensive error handling:

```rust
use savitzky_golay::{SavitzkyGolayError, SavitzkyGolayFilter};

match SavitzkyGolayFilter::new(4, 2) { // Invalid: even window size
    Ok(filter) => { /* use filter */ },
    Err(SavitzkyGolayError::InvalidWindowSize(size)) => {
        println!("Window size {} must be odd", size);
    },
    Err(e) => println!("Error: {}", e),
}
```

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under either of

 * Apache License, Version 2.0, ([LICENSE-APACHE](LICENSE-APACHE) or http://www.apache.org/licenses/LICENSE-2.0)
 * MIT license ([LICENSE-MIT](LICENSE-MIT) or http://opensource.org/licenses/MIT)

at your option.

## References

- Savitzky, A.; Golay, M.J.E. (1964). "Smoothing and Differentiation of Data by Simplified Least Squares Procedures". Analytical Chemistry. 36 (8): 1627–1639.
- Numerical Recipes in C: The Art of Scientific Computing
- [Wikipedia: Savitzky-Golay filter](https://en.wikipedia.org/wiki/Savitzky%E2%80%93Golay_filter)



# Modes
Savitzky-Golay filters require handling signal edges where the moving window extends beyond the data. The modes you mentioned (e.g., from MATLAB's sgolayfilt or SciPy's implementations) specify how to pad or extrapolate the signal at boundaries. These aren't directly part of coefficient computation (which your compute_coefficients function handles) but apply during convolution with the signal.

Here's a brief overview of common modes:

Mirror (or "symmetric"): Reflects the signal at the boundary. For example, if the signal is [a, b, c] and the window needs padding, it becomes [b, a, b, c, b] (mirroring around the edge). This preserves symmetry and is good for smooth signals.

Constant (or "edge"): Pads with the constant value of the nearest edge point. E.g., [a, b, c] becomes [a, a, b, c, c]. Simple but can introduce discontinuities.

Nearest: Repeats the nearest value. Similar to constant but can be more abrupt; often equivalent to constant for edges.

Wrap (or "periodic"): Wraps the signal cyclically. E.g., [a, b, c] becomes [c, a, b, c, a]. Useful for periodic signals but can create artifacts if the signal isn't truly periodic.

Interp (or "linear"): Linearly interpolates or extrapolates values beyond the edges. E.g., extrapolate based on the slope at the boundary. This can reduce edge artifacts but is more computationally intensive.