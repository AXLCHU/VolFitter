#pragma once

#include "Models.hpp"
#include "MarketDataReader.hpp"
#include <vector>
#include <cmath>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include "MathsUtils.hpp"
#include "Fitter.hpp"

// Including NLopt C++ header for COBYLA solver
#include <nlopt.hpp>


void VolFitter::FitSurface(std::vector<VolSlice>& slices, FitterType type) {
    if (slices.empty()) {
        throw std::runtime_error("No market data slices provided for calibration.");
    }

    // 1. Sort slices by expiry in DESCENDING order (longest expiry first)
    std::sort(slices.begin(), slices.end(), [](const VolSlice& a, const VolSlice& b) {
        return a.expiry > b.expiry;
        });

    std::cout << "[VolFitter] Starting slice-by-slice calibration in descending order..." << std::endl;
    std::cout << "[VolFitter] Range: " << slices.front().expiry << " Y down to " << slices.back().expiry << " Y" << std::endl;

    SVIParameters* last_svi = nullptr;
    SSVIParameters* last_ssvi = nullptr;


    // 2. Iterate through each slice and execute the selected optimization
    double penalty_scaler = 40.0;

    for (auto& slice : slices) {
        std::cout << " -> Calibrating slice T = " << slice.expiry << "..." << std::endl;

        if (type == FitterType::SVI) {
            slice.svi_p = FitSingleSliceSVI(slice);
            last_svi = &(slice.svi_p);
            std::cout << "    [SVI Success] a=" << slice.svi_p.a << ", b=" << slice.svi_p.b
                << ", rho=" << slice.svi_p.rho << ", m=" << slice.svi_p.m << ", sigma=" << slice.svi_p.sigma << std::endl;
        }
        else {
            slice.ssvi_p = FitSingleSliceSSVI(slice, last_ssvi, penalty_scaler);
            last_ssvi = &(slice.ssvi_p);
            std::cout << "    [SSVI Success] theta=" << slice.ssvi_p.theta
                << ", phi=" << slice.ssvi_p.phi << ", rho=" << slice.ssvi_p.rho << std::endl;
        }
    }
}


double VolFitter::DeduceAtmVol(const VolSlice& slice) {
    std::vector<std::pair<double, double>> points;
    for (const auto& pt : slice.points) {
        points.push_back({ pt.log_moneyness, pt.implied_vol });
    }
    std::sort(points.begin(), points.end());

    std::vector<double> ks, vols;
    for (const auto& p : points) {
        ks.push_back(p.first);
        vols.push_back(p.second);
    }

    if (ks.size() >= 3) {
        MathUtils::NaturalCubicSpline spline(ks, vols);
        return spline.Evaluate(0.0); // k = 0 signifie ATM
    }

    auto it = std::lower_bound(ks.begin(), ks.end(), 0.0);
    if (it == ks.end()) return vols.back();
    if (it == ks.begin()) return vols.front();
    size_t idx = std::distance(ks.begin(), it);

    return MathUtils::LinearInterpolation(ks[idx - 1], vols[idx - 1], ks[idx], vols[idx], 0.0);
}

// ============================================================================
// SVI CALIBRATION
// ============================================================================

/**
    * NLopt objective function wrapper for SVI Vega-weighted RMSE.
    */
double SviObjective(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    // COBYLA is derivative-free, grad will always be empty
    const VolSlice* slice = static_cast<const VolSlice*>(data);

    SVIParameters p{ x[0], x[1], x[2], x[3], x[4] }; // a, b, rho, m, sigma

    double T = slice->expiry;
    double weighted_error = 0.0;
    double total_weight = 0.0;

    for (const auto& pt : slice->points) {
        double model_w = VolModels::SVI_Raw(pt.log_moneyness, p);
        double model_vol = (model_w > 0.0) ? std::sqrt(model_w / T) : 0.0;
        //            double diff = model_vol - pt.implied_vol;
        double diff = model_w - pt.total_variance;

        double w_vega = pt.vega;
        weighted_error += diff * diff * w_vega;
        total_weight += w_vega;
    }

    return std::sqrt(weighted_error / (total_weight + 1e-9));
}

/**
    * Strict no-arbitrage constraint for SVI Raw parameters.
    * Ensures total variance remains strictly positive everywhere: a + b * sigma * sqrt(1 - rho^2) >= 0
    *          Formatted as g(x) <= 0  ==>  -(a + b * sigma * sqrt(1 - rho^2)) <= 0
    */
double SviArbitrageConstraint(const std::vector<double>& x, std::vector<double>& grad, void* data) {
    double a = x[0];
    double b = x[1];
    double rho = x[2];
    double sigma = x[4];

    double min_variance = a + b * sigma * std::sqrt(1.0 - rho * rho);
    return -min_variance; // Must be <= 0, so min_variance must be >= 0
}

SVIParameters VolFitter::FitSingleSliceSVI(const VolSlice& slice, const SVIParameters* prev_params) {
    // 5 Parameters: a, b, rho, m, sigma
    nlopt::opt opt(nlopt::LN_COBYLA, 5); // LN_COBYLA, LN_BOBYQA

    // Parameter Bounds
    std::vector<double> lb({ 0.0,  0.005, -0.95, -1.0,  0.001 }); // a, b, rho, m, sigma
    std::vector<double> ub({ 0.8,  0.4,  -0.3,  1.0,  0.25 });
    opt.set_lower_bounds(lb);
    opt.set_upper_bounds(ub);

    // Set Objective Function and pass slice data pointer
    opt.set_min_objective(SviObjective, const_cast<VolSlice*>(&slice));

    // Add Non-Linear Arbitrage Constraint (Tolerance set to 1e-8)
    opt.add_inequality_constraint(SviArbitrageConstraint, nullptr, 1e-8);

    // Stopping Criteria
    opt.set_xtol_rel(1e-6);
    opt.set_maxeval(2000);

    // Initial Guess & Step Sizes
    // need to use previous slice as initial guess
    std::vector<double> x0(5); ///

    //opt.set_initial_step(std::vector<double>({ 0.01, 0.02, 0.05, 0.05, 0.02 }));

    if (prev_params != nullptr) {
        // A. On utilise les paramètres de la calibration précédente (maturité plus longue)
        double min_mkt_var = 999.0;
        for (const auto& pt : slice.points) {
            double var = pt.implied_vol * pt.implied_vol * slice.expiry;
            if (var < min_mkt_var) min_mkt_var = var;
        }

        // 2. On adapte x0
        x0[0] = std::max(0.0001, min_mkt_var * 0.9); // 'a' recalibré sur le niveau actuel !
        x0[1] = prev_params->b * 0.8;                // 'b' légèrement réduit pour anticiper des ailes plus plates
        x0[2] = prev_params->rho;                    // 'rho' conservé (le skew reste similaire)
        x0[3] = prev_params->m;                      // 'm' conservé
        x0[4] = std::min(0.10, prev_params->sigma);  // 'sigma' bridé pour éviter un creux trop large
    }
    else {
        // B. C'est la toute première tranche : calcul d'un x0 heuristique basé sur le marché
        double min_mkt_vol = 999.0;
        for (const auto& pt : slice.points) {
            if (pt.implied_vol < min_mkt_vol) min_mkt_vol = pt.implied_vol;
        }
        double estimated_min_var = min_mkt_vol * min_mkt_vol * slice.expiry;

        ////
        x0[0] = std::max(0.001, estimated_min_var * 0.8); // a : un peu en dessous de la variance min
        x0[1] = 0.15;                                     // b : pente standard
        x0[2] = -0.70;                                    // rho : skew négatif marqué (SP500)
        x0[3] = 0.0;                                      // m : centré sur le forward
        x0[4] = 0.05;                                     // sigma : courbure initiale
    }

    // Sécurisation de x0 par rapport aux bornes au cas où prev_params flirterait avec les limites
    for (int i = 0; i < 5; ++i) {
        if (x0[i] < lb[i]) x0[i] = lb[i];
        if (x0[i] > ub[i]) x0[i] = ub[i];
    }

    double minf;
    try {
        nlopt::result result = opt.optimize(x0, minf);
    }
    catch (const std::exception& e) {
        std::cerr << " [SVI Warning] Optimization exception on T=" << slice.expiry << " : " << e.what() << std::endl;
    }

    SVIParameters out_p;
    out_p.a = x0[0]; out_p.b = x0[1]; out_p.rho = x0[2]; out_p.m = x0[3]; out_p.sigma = x0[4];
    return out_p;
}


// ============================================================================
// SSVI CALIBRATION
// ============================================================================

struct SsviArgs {
    const VolSlice* slice;
    FitterType type;
    double theta;
    const SSVIParameters* next_slice_params;
    double penalty_scaler;
};

double SsviVolPenaltyTerm(double var, double t, double next_var, double s) {
    if (t <= 0) return 0.0;
    double diff = (var * t - next_var) / t;
    return s * std::max(diff, 0.0);
}

double SsviVolErrors(const std::vector<double>& params, std::vector<double>& grad, void* data) {
    
    SsviArgs* args = static_cast<SsviArgs*>(data);
    const VolSlice* slice = args->slice;

    double theta = args->theta;
    double T = args->slice->expiry;
    double atm = std::sqrt(theta / T);
    SSVIParameters p{ args->theta, params[1], params[0] };
    double f = 0.0;

    for (size_t i = 0; i < slice->points.size(); ++i) {
        const auto& pt = slice->points[i];
        double var_model = VolModels::SSVI(pt.log_moneyness, p);
        double iv = pt.implied_vol;
        double diff_variance = var_model - (iv * iv);
        double error_term = (diff_variance * diff_variance) * pt.vega;
        f += error_term;

        if (args->next_slice_params != nullptr) {
            double theta_next = args->next_slice_params->theta;
            double rho_next = args->next_slice_params->rho;
            double phi_next = args->next_slice_params->phi;

            double tvar_next_slice_i = VolModels::SSVI(pt.log_moneyness, SSVIParameters{ theta_next, phi_next, rho_next });

            f += SsviVolPenaltyTerm(var_model, slice->expiry, tvar_next_slice_i, args->penalty_scaler) * pt.vega;
        }
    }

    return f;
}

double SsviButterflyConstraint(const std::vector<double>& params, std::vector<double>& grad, void* data) {
    
    double s = params[0]; // rho
    double c = params[1]; // phi

    double v = *static_cast<double*>(data);

    double c0_temp1 = (std::pow(1.0 - ((1.0 / 8.0) * v), 2)) + v;
    double c0_temp2 = (5.0 - ((1.0 / 8.0) * v)) / c0_temp1;
    double c0 = c0_temp2 + std::sqrt(std::pow(c0_temp2, 2) - (1.0 / c0_temp1));

    double s0 = std::sqrt(4.0 - v);

    return 1.0 - std::abs(s) / s0 - c / c0;
}

SSVIParameters VolFitter::FitSingleSliceSSVI(const VolSlice& slice, const SSVIParameters* prev_params, double penalty_scaler) {

    double atm_vol = VolFitter::DeduceAtmVol(slice);
    double theta = atm_vol * atm_vol * slice.expiry;

    SsviArgs args;
    args.slice = &slice;
    args.theta = theta;
    args.penalty_scaler = penalty_scaler;

    SSVIParameters next_params;
    if (prev_params != nullptr) {
        next_params.theta = prev_params->theta;
        next_params.rho = prev_params->rho;
        next_params.phi = prev_params->phi;
        args.next_slice_params = &next_params;
    }
    else {
        args.next_slice_params = nullptr;
    }

    std::vector<double> xn(2); // Last success guess (xn)
    if (prev_params != nullptr) {
        xn[0] = prev_params->rho;
        xn[1] = prev_params->phi;
    }
    else {
        xn[0] = -0.5; // S2_GUESS
        xn[1] = 0.5;  // C2_GUESS
    }

    std::vector<double> x0 = { -0.5, 0.5 }; // Initial guess par défaut (x0)

    std::vector<nlopt::algorithm> optimizers = { nlopt::LD_SLSQP, nlopt::LN_COBYLA };
    std::vector<std::vector<double>> guesses = { xn, x0 };

    std::vector<double> lb = { -0.999, 0.05 }; // s2_min (rho), c2_min (phi)
    std::vector<double> ub = { 0.999, 10.0 };

    std::vector<double> best_x = x0;
    bool success = false;

    for (const auto& guess : guesses) {
        for (const auto& algo : optimizers) {
            nlopt::opt opt(algo, 2);
            opt.set_lower_bounds(lb);
            opt.set_upper_bounds(ub);
            opt.set_min_objective(SsviVolErrors, &args);

            double v_param = theta; // atm^2 * t
            opt.add_inequality_constraint(SsviButterflyConstraint, &v_param, 1e-8);

            opt.set_xtol_rel(1e-12);
            opt.set_maxeval(1500);

            std::vector<double> x = guess;
            double minf;
            try {
                if (opt.optimize(x, minf) == nlopt::SUCCESS) {
                    best_x = x;
                    success = true;
                    break;
                }
            }
            catch (...) {
            }
        }
        if (success) break;
    }

    SSVIParameters out_p;
    out_p.theta = theta;
    out_p.rho = best_x[0];
    out_p.phi = best_x[1];
    return out_p;
}
