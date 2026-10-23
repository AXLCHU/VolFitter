#include "ArbitrageChecker.hpp"
#include "Models.hpp"
#include <cmath>
#include <algorithm>
#include <iostream>


std::vector<ArbitrageViolation> ArbitrageChecker::CheckMarketData(const std::vector<VolSlice>& slices) {
    std::vector<ArbitrageViolation> violations;
    if (slices.empty()) return violations;

    // ------------------------------------------------------------------------
    // BUTTERFLY ARB
    // ------------------------------------------------------------------------
    for (const auto& slice : slices) {
        if (slice.points.size() < 3) continue;

		// Sort by increasing strike to ensure correct order for finite difference
        auto sorted_points = slice.points;
        std::sort(sorted_points.begin(), sorted_points.end(), [](const MarketPoint& a, const MarketPoint& b) {
            return a.strike < b.strike;
            });

        for (size_t i = 1; i < sorted_points.size() - 1; ++i) {
            double K_prev = sorted_points[i - 1].strike;
            double K_curr = sorted_points[i].strike;
            double K_next = sorted_points[i + 1].strike;

            double w_prev = sorted_points[i - 1].total_variance;
            double w_curr = sorted_points[i].total_variance;
            double w_next = sorted_points[i + 1].total_variance;

            double h1 = K_curr - K_prev;
            double h2 = K_next - K_curr;
            if (h1 <= 0.0 || h2 <= 0.0) continue;

            double w_prime_prime_K = 2.0 * ((w_next / h2) - w_curr * (1.0 / h1 + 1.0 / h2) + (w_prev / h1)) / (h1 + h2);

            if (w_prime_prime_K < -1e-7) {
                violations.push_back({ "Butterfly", slice.expiry, sorted_points[i].strike, w_prime_prime_K });
            }
        }
    }

    // ------------------------------------------------------------------------
    // CALENDAR ARB
    // ------------------------------------------------------------------------
    auto chrono_slices = slices;
    std::sort(chrono_slices.begin(), chrono_slices.end(), [](const VolSlice& a, const VolSlice& b) {
        return a.expiry < b.expiry;
        });

    for (size_t t = 0; t < chrono_slices.size() - 1; ++t) {
        const auto& slice_short = chrono_slices[t];
        const auto& slice_long = chrono_slices[t + 1];

        for (const auto& p_short : slice_short.points) {
            for (const auto& p_long : slice_long.points) {
                if (std::abs(p_short.strike - p_long.strike) < 1e-2) {

                    if (p_long.total_variance < p_short.total_variance - 1e-6) {
                        double diff = p_long.total_variance - p_short.total_variance;
                        violations.push_back({ "Calendar", slice_long.expiry, p_short.strike, diff });
                    }
                    break;
                }
            }
        }
    }

    return violations;
}



std::vector<ArbitrageViolation> ArbitrageChecker::CheckSurface(const std::vector<VolSlice>& slices, FitterType type) {
    
    std::vector<ArbitrageViolation> violations;
    if (slices.empty()) return violations;

    // ------------------------------------------------------------------------
    // BUTTERFLY ARBITRAGE CHECK (Durrleman / Gatheral Density Condition)
    // ------------------------------------------------------------------------
    // The condition g(k) >= 0 ensures non-negative local density, where:
    // g(k) = (1 - (k * w') / (2 * w))^2 - (w' / 4) * (1 / w + 1 / 4) + w'' / 2

    // Define a fine log-moneyness grid to catch subtle wing flips
    std::vector<double> k_grid;
    for (double k = -3.0; k <= 3.0; k += 0.01) {
        k_grid.push_back(k);
    }
    double dk = 0.0001; // Step for central finite-difference derivatives

    for (const auto& slice : slices) {
        for (double k : k_grid) {
            // Evaluate continuous model values at k, k-dk, k+dk using the calibrated parameters
            double w = 0.0, w_plus = 0.0, w_minus = 0.0;

            if (type == FitterType::SVI) {
                w = VolModels::SVI_Raw(k, slice.svi_p);
                w_plus = VolModels::SVI_Raw(k + dk, slice.svi_p);
                w_minus = VolModels::SVI_Raw(k - dk, slice.svi_p);
            }
            else {
                w = VolModels::SSVI(k, slice.ssvi_p);
                w_plus = VolModels::SSVI(k + dk, slice.ssvi_p);
                w_minus = VolModels::SSVI(k - dk, slice.ssvi_p);
            }

            // Numerical derivatives with respect to log-moneyness k
            double w_prime = (w_plus - w_minus) / (2.0 * dk);
            double w_prime_prime = (w_plus - 2.0 * w + w_minus) / (dk * dk);

            // Compute Gatheral's g(k) metric
            double term1 = 1.0 - (k * w_prime) / (2.0 * w);
            double term2 = (w_prime * w_prime / 4.0) * (1.0 / w + 0.25);
            double g_k = (term1 * term1) - term2 + (w_prime_prime / 2.0);

            // If g(k) < 0, butterfly arbitrage exists (negative probability density)
            if (g_k < -1e-6) { // Small numerical tolerance threshold
                // Reconstruct approximate absolute strike: K = F * exp(k)
                // Assuming a generic reference forward of 100 for reporting context
                double estimated_strike = 100.0 * std::exp(k);
                violations.push_back({ "Butterfly", slice.expiry, estimated_strike, g_k });
            }
        }
    }

    // ------------------------------------------------------------------------
    // CALENDAR ARBITRAGE CHECK (Time Monotonicity)
    // ------------------------------------------------------------------------
    // Calendar arbitrage-free condition requires: w(k, T_2) >= w(k, T_1) for T_2 > T_1
    // Total variance must be a non-decreasing function of time at any fixed k.

    // Clone and sort slices by chronological CHRONOLOGICAL order (shortest maturity first)
    auto chronological_slices = slices;
    std::sort(chronological_slices.begin(), chronological_slices.end(), [](const VolSlice& a, const VolSlice& b) {
        return a.expiry < b.expiry;
        });

    for (size_t i = 0; i < chronological_slices.size() - 1; ++i) {
        const auto& current_slice = chronological_slices[i];
        const auto& next_slice = chronological_slices[i + 1];

        for (double k : k_grid) {
            double w1 = 0.0;
            double w2 = 0.0;

            if (type == FitterType::SVI) {
                w1 = VolModels::SVI_Raw(k, current_slice.svi_p);
                w2 = VolModels::SVI_Raw(k, next_slice.svi_p);
            }
            else {
                w1 = VolModels::SSVI(k, current_slice.ssvi_p);
                w2 = VolModels::SSVI(k, next_slice.ssvi_p);
            }

            // Total variance dropped over a positive time interval
            if (w2 < w1 - 1e-6) {
                double estimated_strike = 100.0 * std::exp(k);
                double delta_w = w2 - w1;
                violations.push_back({ "Calendar", next_slice.expiry, estimated_strike, delta_w });
            }
        }
    }

    return violations;
}


void UpdatePointWithModel(MarketPoint& pt, double F, double T, const VolSlice& slice, FitterType type) {

    double k = std::log(pt.strike / F);
    pt.log_moneyness = k;

    double total_var = 0.0;
    if (type == FitterType::SVI) {
        total_var = VolModels::SVI_Raw(k, slice.svi_p);
    }
    else {
        total_var = VolModels::SSVI(k, slice.ssvi_p);
    }

    if (total_var <= 0.0) {
        pt.total_variance = 0.0;
        pt.implied_vol = 0.0;
    }
    else {
        pt.total_variance = total_var;
        pt.implied_vol = std::sqrt(total_var / T); // Vol = sqrt(w / T)
    }
}