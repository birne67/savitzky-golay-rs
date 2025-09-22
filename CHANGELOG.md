# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.1.0] - 2025-09-22

### Added
- Initial implementation of Savitzky-Golay filter
- Core mathematical algorithms for coefficient calculation using least squares
- Support for different window sizes and polynomial orders
- Multiple boundary handling strategies (Constant, Mirror, Wrap, Zero, Shrink)
- Numerical differentiation capabilities (any order)
- Comprehensive error handling with custom error types
- Performance optimizations with coefficient caching
- Extensive test suite with mathematical validation
- Documentation with examples and mathematical background
- Convenience functions for common use cases (`smooth()`, `derivative()`)
- Multiple examples demonstrating usage and performance

### Features
- High-performance filtering with ~3ms for 10,000 points
- Polynomial preservation up to the specified order
- Zero-copy operations where possible
- Thread-safe operations
- Comprehensive API documentation
- Support for f64 data types

[Unreleased]: https://github.com/yourusername/savitzky-golay/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/yourusername/savitzky-golay/releases/tag/v0.1.0