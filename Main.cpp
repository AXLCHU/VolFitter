#include "Models.hpp"
#include "MarketDataReader.hpp"
#include "Fitter.hpp"
#include "ArbitrageChecker.hpp"
#include <iostream>
#include <vector>
#include <string>
#include "VolSurfaceEvaluator.hpp"


int main() {

    // ------------------------------------------------------------------------
    // LOADING MARKET DATA
    // ------------------------------------------------------------------------
    std::string forward_csv = "C:\\Users\\axelc\\OneDrive\\Documents\\github\\VolFitter\\inputs\\20250630\\forward_sp500.csv";
    std::string vol_csv = "C:\\Users\\axelc\\OneDrive\\Documents\\github\\VolFitter\\inputs\\20250630\\vols_sp500.csv";
    std::string output_path = "C:\\Users\\axelc\\OneDrive\\Documents\\github\\VolFitter\\outputs\\20250630\\";

    std::cout << "\n[Step 1] Loading Market Data..." << std::endl;

    std::cout << "=== Loading Forwards ===" << std::endl;
    std::map<int, double> forward_curve;
    try {
        forward_curve = MarketDataReader::ReadForwardCurve(forward_csv);
        std::cout << " -> Success: " << forward_curve.size() << " Forwards maturities loaded.\n";
    }
    catch (const std::exception& e) {
        std::cerr << " -> Error: " << e.what() << std::endl;
        return 1;
    }

    std::cout << "\n=== Loading Volatility Surface ===" << std::endl;
    std::vector<VolSlice> surface_data;
    try {
        surface_data = MarketDataReader::ReadCSV(vol_csv, forward_curve);
        std::cout << " -> Success.\n";
    }
    catch (const std::exception& e) {
        std::cerr << " -> Error: " << e.what() << std::endl;
        return 1;
    }

    // Check arb
    std::cout << "\n[Step 1.5] Check Market Data Arb..." << std::endl;
    std::vector<ArbitrageViolation> market_violations = ArbitrageChecker::CheckMarketData(surface_data);

    if (market_violations.empty()) {
        std::cout << " -> No Arb.\n";
    }
    else {
        std::cout << " -> ATTENTION : " << market_violations.size() << " arb in input:\n";
        for (const auto& v : market_violations) {
            std::cout << "    * [" << v.type << "] at T = "
                << v.expiry << "Y | Strike = " << v.strike
                << " | Value = " << v.metric_val << std::endl;
        }
    }


    // ------------------------------------------------------------------------
    // SELECTING MODEL & RUNNING OPTIMIZATION
    // ------------------------------------------------------------------------

    //FitterType selected_model = FitterType::SVI;
    //FitterType selected_model = FitterType::SSVI;
    //FitterType selected_model = FitterType::SVI_JW;
    FitterType selected_model = FitterType::SVI_wings;

    switch (selected_model) {
    case FitterType::SVI:
        std::cout << "\n[Step 2] Executing SVI Volatility Fitter..." << std::endl;
        break;
    case FitterType::SSVI:
        std::cout << "\n[Step 2] Executing SSVI Volatility Fitter..." << std::endl;
        break;
    case FitterType::SVI_JW:
        std::cout << "\n[Step 2] Executing SVI Jump-Wings (SVI-JW) Volatility Fitter..." << std::endl;
        break;
    case FitterType::SVI_wings:
        std::cout << "\n[Step 2] Executing SVI Power-Law Wings Volatility Fitter..." << std::endl;
        break;
    }

    try {
        VolFitter::FitSurface(surface_data, selected_model);
        std::cout << "\n -> Optimization process completed successfully." << std::endl;
    }
    catch (const std::exception& e) {
        std::cerr << "\n -> ERROR: Optimization failed. " << e.what() << std::endl;
        return 1;
    }

    // ------------------------------------------------------------------------
    // POST-CALIBRATION ARBITRAGE DIAGNOSTIC
    // ------------------------------------------------------------------------
    std::cout << "\n[Step 3] Running surface arbitrage check..." << std::endl;

    // Check at the input data points only
    std::vector<VolSlice> fitted_surface_data = VolSurfaceEvaluator::GenerateFittedSurface(surface_data, forward_curve, selected_model);

    // Output fitted vols
    std::string output_csv = output_path + (selected_model == FitterType::SVI ? "svi_vols.csv" : 
                                            selected_model == FitterType::SSVI ? "ssvi_vols.csv" : 
                                            selected_model == FitterType::SVI_JW ? "svijw_vols.csv" : "svi_wings_vols.csv");

    VolSurfaceEvaluator::ExportToCSV(output_csv, fitted_surface_data);

    // On lance le checker sur l'objet converti
    std::vector<ArbitrageViolation> violations = ArbitrageChecker::CheckMarketData(fitted_surface_data);

    // Scan the generated surface over a dense grid to locate butterfly or calendar violations
    //std::vector<ArbitrageViolation> violations = ArbitrageChecker::CheckSurface(surface_data, selected_model);

    if (violations.empty()) {
        std::cout << " -> EXCELLENT: The fitted surface is 100% Static Arbitrage Free!" << std::endl;
    }
    else {
        std::cout << " -> WARNING: Detected " << violations.size() << " arbitrage violation(s) on the surface:" << std::endl;

        for (const auto& v : violations) {
            std::cout << "    * [" << v.type << " Arbitrage] found at Expiry T = "
                << v.expiry << " Y | Strike K = " << v.strike << std::endl;
                //<< " | Metric Value = " << v.metric_val << std::endl;
        }
    }

    std::cout << "\nEnd of Volatility Surface calibration." << std::endl;

    return 0;
}
