#pragma once
#include <vector>
#include <string>
#include <map>
#include "Models.hpp"


struct MarketPoint {
    double strike;
    double implied_vol;
    double vega; // Calculated when reading using Black-Scholes
    double log_moneyness;
    double total_variance;
};

struct VolSlice {
    double expiry;
    int expiry_days;
    std::vector<MarketPoint> points;
    SVIParameters svi_p;
    SSVIParameters ssvi_p;
    SVIJWParameters svijw_p;
	PowerLawWingsParameters svi_wings_p;
};

class MarketDataReader {
public:
    static std::map<int, double> ReadForwardCurve(const std::string& filepath);
    // Read CSV and convert strikes/IV to log-moneyness/total_variance
    static std::vector<VolSlice> ReadCSV(const std::string& filepath, const std::map<int, double>& forward_curve);
};
