#pragma once
#include <cmath>

enum class FitterType { SVI, SSVI };

struct SVIParameters { double a, b, rho, m, sigma; };
struct SSVIParameters { double theta, phi, rho; };

namespace VolModels {

    inline double SVI_Raw(double k, const SVIParameters& p) {
        return p.a + p.b * (p.rho * (k - p.m) + std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma));
    }

/*    inline double SSVI(double k, const SSVIParameters& p) {
        return 0.5 * p.theta * (1.0 + p.rho * p.phi * k +
            std::sqrt((p.phi * k + p.rho) * (p.phi * k + p.rho) + (1.0 - p.rho * p.rho)));
    }*/

    inline double SSVI(double k, const SSVIParameters& p) {
        
        double total_atm_vol = std::sqrt(p.theta);

        if (total_atm_vol <= 0.0) return 0.0;

        double z = k / total_atm_vol;

        double term1 = 0.5 * (1.0 + p.rho * z);
        double inner_sqrt = 0.25 * std::pow(1.0 + p.rho * z, 2) + 0.5 * p.phi * std::pow(z, 2);

        if (inner_sqrt < 0.0) inner_sqrt = 0.0;

        return p.theta * (term1 + std::sqrt(inner_sqrt));
    }
}
