#pragma once
#include <cmath>
#include <algorithm>

enum class FitterType { SVI, SSVI, SVI_JW, SVI_wings };

struct SVIParameters { double a, b, rho, m, sigma; };
struct SVIJWParameters { double vt, psi, p, c, vt_til; };
struct SSVIParameters { double theta, phi, rho; };
struct IVPParameters { double a, b, rho, m, sigma, gamma; }; // IVP model (SVI with global power)

// Power-Law Wing SVI
struct PowerLawWingsParameters {
    SVIParameters svi_core; 
    double k_L, A_L, gamma_L; // Put 
    double k_R, A_R, gamma_R; // Call 
};


namespace VolModels {

    inline double SVI_Raw(double k, const SVIParameters& p) {
        return p.a + p.b * (p.rho * (k - p.m) + std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma));
    }

    inline double SVI_Raw_Derivative(double k, const SVIParameters& p) {
        double sqrt_term = std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma);
        if (sqrt_term < 1e-9) sqrt_term = 1e-9;
        return p.b * (p.rho + (k - p.m) / sqrt_term);
    }

/*    inline double SSVI(double k, const SSVIParameters& p) {
        return 0.5 * p.theta * (1.0 + p.rho * p.phi * k +
            std::sqrt((p.phi * k + p.rho) * (p.phi * k + p.rho) + (1.0 - p.rho * p.rho)));
    }*/


    // Convert Jump-Wings in Raw SVI parameters 
    inline SVIParameters ConvertJwToRaw(const SVIJWParameters& jw, double t) {
        double wt = jw.vt * t;
        if (wt < 1e-9) wt = 1e-9;

        double beta = jw.p + jw.c;
        double alpha = jw.p - jw.c;

        SVIParameters raw;

        raw.b = 0.5 * std::sqrt(wt) * beta;
        raw.rho = (beta != 0.0) ? (alpha / beta) : 0.0;

        if (raw.rho >= 1.0)  raw.rho = 0.999;
        if (raw.rho <= -1.0) raw.rho = -0.999;

        if (raw.b * (1.0 - raw.rho * raw.rho) == 0.0) {
            raw.m = 0.0;
            raw.sigma = 0.05;
        }
        else {
            double denom = raw.b * std::sqrt(1.0 - raw.rho * raw.rho);
            raw.m = (jw.vt - jw.vt_til) / (denom * t);

            double g_atm = (2.0 * jw.psi * std::sqrt(wt)) / raw.b;
            double sqrt_term = 1.0 - raw.rho * raw.rho;
            raw.sigma = std::abs(raw.m) * std::sqrt(std::max(1.0 / (g_atm * g_atm - raw.rho * raw.rho) - 1.0, 1e-6));
        }

        raw.a = jw.vt_til * t - raw.b * raw.sigma * std::sqrt(1.0 - raw.rho * raw.rho);

        return raw;
    }

    // SVI-JW tvar at k and t
    inline double SVI_JW(double k, const SVIJWParameters& p, double t) {
        SVIParameters raw = ConvertJwToRaw(p, t);
        return SVI_Raw(k, raw);
    }

    inline double SSVI(double k, const SSVIParameters& p) {
        
        double total_atm_vol = std::sqrt(p.theta);

        if (total_atm_vol <= 0.0) return 0.0;

        double z = k / total_atm_vol;

        double term1 = 0.5 * (1.0 + p.rho * z);
        double inner_sqrt = 0.25 * std::pow(1.0 + p.rho * z, 2) + 0.5 * p.phi * std::pow(z, 2);

        if (inner_sqrt < 0.0) inner_sqrt = 0.0;

        return p.theta * (term1 + std::sqrt(inner_sqrt));
    }

    inline double IVP_TotalVariance(double k, const IVPParameters& p) {
        double inner = p.rho * (k - p.m) + std::sqrt((k - p.m) * (k - p.m) + p.sigma * p.sigma);
        if (inner < 0.0) inner = 0.0; 
        return p.a + p.b * std::pow(inner, 2.0 * p.gamma);
    }

    inline double SVI_With_PowerLawWings(double k, const PowerLawWingsParameters& p) {
        
        if (k < p.k_L) {
            return p.A_L * std::pow(std::abs(k), p.gamma_L);
        }
        else if (k > p.k_R) {
            return p.A_R * std::pow(std::abs(k), p.gamma_R);
        }
        else {
            return SVI_Raw(k, p.svi_core);
        }
        //return SVI_Raw(k, p.svi_core);
    }


}
