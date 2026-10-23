#pragma once

#include "MarketDataReader.hpp"
#include "Models.hpp"
#include <vector>
#include <string>


struct ArbitrageViolation {
    std::string type;   // "Butterfly" or "Calendar"
    double expiry;      // Expiry time T where the violation occurs
    double strike;      // Strike price K where the violation occurs
    double metric_val;  // The underlying math value that failed (e.g., g(k) < 0 or dw < 0)
};

class ArbitrageChecker {
    public:
        static std::vector<ArbitrageViolation> CheckMarketData(const std::vector<VolSlice>& slices);
        static std::vector<ArbitrageViolation> CheckSurface(const std::vector<VolSlice>& slices, FitterType type);
};
