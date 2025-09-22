use savitzky_golay::compute_coefficients;

fn main() {
    // Diagnostic: print derivative coefficients for window=7, poly_order=4, derivative=1
    let window = 7;
    let poly = 4;
    let deriv = 1;

    let coeffs = compute_coefficients(window, poly, deriv).expect("compute coeffs failed");
    println!("coeffs (window={}, poly={}, deriv={}):", window, poly, deriv);
    for (i, c) in coeffs.iter().enumerate() {
        println!("  [{}] = {:.12}", i, c);
    }

    // Create test signal f(x)=x^4 sampled at step=0.1
    let step = 0.1;
    let n = 20;
    let data: Vec<f64> = (0..n).map(|i| {
        let x = i as f64 * step;
        x.powi(4)
    }).collect();

    // Apply convolution at a chosen interior index, e.g., i=10
    let center = 10usize;
    let half = window / 2;
    let mut conv = 0.0;
    for j in 0..window {
        let idx = (center + j).saturating_sub(half);
        conv += coeffs[j] * data[idx];
    }

    println!("\nAt center index {}: conv = {:.12}", center, conv);

    // Expected derivative at x=center*step: f'(x)=4 x^3 / step
    let x = center as f64 * step;
    let expected = 4.0 * x.powi(3) / step;
    println!("expected derivative (4 x^3 / step) = {:.12}", expected);
}
