
# Equity Volatility Surface Calibration Engine in C++

Implementation and comparison of SVI, SSVI and Power-Law Wings volatility models, calibrated on SPX volatility data.

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

- **SVI:** Parametric IV smile with 5 parameters, fit slice by slice, good fit but prone to arbitrage.
- **SSVI:** 2-3 Parameters Surface parameterisation across maturities and arbitrage-free constraints.
- **Power-Law Wings:** 6 Parameters Concave curvature on extrem wings, ensuring compliance with Roger Lee's moment formula.

## #Arbitrage-Removal-Procedure

To ensure a valid risk-neutral density, our pipeline enforces strict static arbitrage filters:

1. **Calendar Spread Filter:** Evaluates total variance monotonicity \(\partial_t w(k, t) \geq 0\). Violations are resolved by shifting local variances upward or enforcing the Gatheral-Jacquier monotonic mapping condition during multi-slice calibration.
2. **Butterfly Arbitrage (Density) Correction:** Computes the implied risk-neutral density g(k). If a local dip creates g(k) < 0:
   * **SSVI:** Analytically bounded using θ φ(θ)(1+|ρ|) < 4.
   * **SVI & ORC:** Solved via sequential quadratic programming (SQP) that penalizes the negative regions of the density function until g(k) ≥ 0 across the defined moneyness grid.

## #Results

Models are compared based on Calibration error, Convergence and parameter stability and Static-arbitrage compliance

* **Fitting Quality:** **ORC Wing** achieves the lowest overall RMSE on extreme OTM wings due to its dedicated curvature parameters, outperforming Raw SVI.
* **Robustness:** **SSVI** provides the most stable surface across time-steps and guarantees zero arbitrage violations natively, though it slightly underfits highly asymmetric skews.
* **Computational Cost:** Calibration speed is highest for SSVI (lowest parameter count), while ORC Wing requires 1.5x more CPU cycles to converge under the arbitrage removal penalization layer.

## Future Work

- Spline and Radial Basis Function (RBF) surface representations