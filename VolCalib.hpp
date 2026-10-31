#pragma once
#include "Models.hpp"


namespace VolModels {

    inline double SVI_Raw_SecondDerivative(double k, const SVIParameters& p) {
        double diff = k - p.m;
        double square_term = diff * diff + p.sigma * p.sigma;
        double sqrt_term = std::sqrt(square_term);
        if (sqrt_term < 1e-9) return 0.0;
        return p.b * (p.sigma * p.sigma) / (square_term * sqrt_term);
    }

    inline double Wing_Objective_Function(double k, double gamma, const SVIParameters& svi) {
        return k * SVI_Raw_Derivative(k, svi) - gamma * SVI_Raw(k, svi);
    }

    inline double Wing_Objective_Derivative(double k, double gamma, const SVIParameters& svi) {
        return (1.0 - gamma) * SVI_Raw_Derivative(k, svi) + k * SVI_Raw_SecondDerivative(k, svi);
    }

    inline double Find_Wing_Threshold(double k_init, double gamma, const SVIParameters& svi,
        int max_iter = 100, double tol = 1e-7) {
        double k = k_init;
        for (int i = 0; i < max_iter; ++i) {
            double f = Wing_Objective_Function(k, gamma, svi);
            double df = Wing_Objective_Derivative(k, gamma, svi);
            if (std::abs(df) < 1e-9) break;
            double delta = f / df;
            k -= delta;
            if (std::abs(delta) < tol) return k;
        }
        return k;
    }

    inline double Compute_Wing_Amplitude(double k_seuil, double gamma, const SVIParameters& svi) {
        double abs_k = std::abs(k_seuil);

        if (abs_k < 1e-3) {
            abs_k = 1e-3;
        }
        return SVI_Raw(k_seuil, svi) / std::pow(abs_k, gamma);
    }

    /**
     * @brief Point d'entrée principal pour calibrer la structure à ailes concaves
     */
    inline PowerLawWingsParameters Calibrate_PowerLaw_Wings(const SVIParameters& svi,
        double gamma_L, double gamma_R) {
        
        PowerLawWingsParameters p;
        p.svi_core = svi;
        p.gamma_L = gamma_L;
        p.gamma_R = gamma_R;

        double m = svi.m;
        double shift = std::max(svi.sigma, 0.1);

        // Points de départ de la recherche (Négatif pour le Put, Positif pour le Call)
        double k_L_init = m - 0.2 * shift;
        double k_R_init = m + 0.2 * shift;

        p.k_L = Find_Wing_Threshold(k_L_init, gamma_L, svi);
        p.k_R = Find_Wing_Threshold(k_R_init, gamma_R, svi);

        p.A_L = Compute_Wing_Amplitude(p.k_L, gamma_L, svi);
        p.A_R = Compute_Wing_Amplitude(p.k_R, gamma_R, svi);

        return p;
    }
}
