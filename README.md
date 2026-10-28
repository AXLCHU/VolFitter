
# Equity Volatility Surface Calibration Engine (C++)

C++ implementation and comparison of SVI, SSVI and SVI-JW volatility models, calibrated on SPX vol data.

## Features

- Nonlinear parameter calibration and optimisation
- Sequential calibration across maturities
- Static-arbitrage constraints and diagnostics
- Comparison of calibration accuracy and stability
- Performance benchmarking

## Calibration Workflow

1. Prepare SPX implied volatility data.
2. Initialise model parameters.
3. Calibrate using constrained nonlinear optimisation.
4. Validate butterfly and calendar arbitrage conditions.
5. Compare model fit, convergence and runtime.

## Models

- **SVI:** Parametric implied volatility smile.
- **SSVI:** Surface parameterisation across maturities.
- **SVI-JW:** SVI parameterisation using interpretable smile characteristics.

## Benchmarks

Evaluate each model based on:

- Calibration error
- Convergence and parameter stability
- Runtime
- Static-arbitrage compliance

## Future Work

- Spline and radial basis function (RBF) surface representations
- Robust calibration under noisy or sparse market data