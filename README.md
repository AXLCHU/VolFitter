
# Equity Volatility Surface Calibration Engine (C++)

C++ implementation and comparison of SVI, SSVI and SVI-JW volatility models, calibrated on SPX vol data.

## Features

- Nonlinear parameter calibration and optimisation
- Sequential calibration across maturities
- Static-arbitrage constraints and diagnostics

## Calibration Workflow

1. Initialise model parameters.
2. Calibrate using constrained nonlinear optimisation.
3. Validate butterfly and calendar arbitrage conditions.
4. Compare model fit, convergence and runtime.

## Model testing

- **SVI:** Parametric IV smile.
- **SSVI:** Surface parameterisation across maturities and arbitrage-free constraints.
- **SVI-JW:** SVI parameterisation with wings characteristics.

## Benchmarks

Evaluate each model based on:

- Calibration error
- Convergence and parameter stability
- Static-arbitrage compliance

## Future Work

- Spline and radial basis function (RBF) surface representations
