#pragma once
#include <vector>
#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace MathUtils {

    // Classe encapsulée pour la Spline Cubique Naturelle
    class NaturalCubicSpline {
    private:
        std::vector<double> x, y, b, c, d;
    public:
        NaturalCubicSpline() {}
        NaturalCubicSpline(const std::vector<double>& x_in, const std::vector<double>& y_in) {
            x = x_in; y = y_in;
            size_t n = x.size();
            if (n < 3) throw std::runtime_error("Spline nécessite au moins 3 points.");

            std::vector<double> h(n - 1), alpha(n - 1, 0.0);
            for (size_t i = 0; i < n - 1; ++i) h[i] = x[i + 1] - x[i];

            std::vector<double> l(n, 1.0), mu(n, 0.0), z(n, 0.0);
            c.resize(n, 0.0); b.resize(n - 1); d.resize(n - 1);

            for (size_t i = 1; i < n - 1; ++i) {
                alpha[i] = (3.0 / h[i]) * (y[i + 1] - y[i]) - (3.0 / h[i - 1]) * (y[i] - y[i - 1]);
                l[i] = 2.0 * (x[i + 1] - x[i - 1]) - h[i - 1] * mu[i - 1];
                mu[i] = h[i] / l[i];
                z[i] = (alpha[i] - h[i - 1] * z[i - 1]) / l[i];
            }

            for (int i = static_cast<int>(n) - 2; i >= 0; --i) {
                c[i] = z[i] - mu[i] * c[i + 1];
                b[i] = (y[i + 1] - y[i]) / h[i] - h[i] * (c[i + 1] + 2.0 * c[i]) / 3.0;
                d[i] = (c[i + 1] - c[i]) / (3.0 / h[i]);
            }
        }

        double Evaluate(double val) const {
            if (x.empty()) return 0.0;
            auto it = std::lower_bound(x.begin(), x.end(), val);
            size_t idx = std::distance(x.begin(), it);
            if (idx == 0) idx = 1;
            if (idx >= x.size()) idx = x.size() - 1;
            size_t i = idx - 1;
            double dx = val - x[i];
            return y[i] + b[i] * dx + c[i] * dx * dx + d[i] * dx * dx * dx;
        }
    };

    // Fonction d'interpolation linéaire classique
    inline double LinearInterpolation(double x1, double y1, double x2, double y2, double x_target) {
        if (std::abs(x2 - x1) < 1e-9) return y1;
        return y1 + (x_target - x1) * (y2 - y1) / (x2 - x1);
    }
}
