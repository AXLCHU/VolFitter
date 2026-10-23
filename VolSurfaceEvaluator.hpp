#pragma once
#include "Models.hpp"
#include "ArbitrageChecker.hpp"
#include <vector>
#include <map>
#include <string>
#include <fstream>
#include <iomanip>
#include <cmath>
#include <iostream>
#include <set>
#include <algorithm>


class VolSurfaceEvaluator {

    private:
        // Calculate total variance and vol for a point 
        static void UpdatePointWithModel(MarketPoint& pt, double T, const VolSlice& slice, FitterType type) {
            double k = pt.log_moneyness;

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
                pt.implied_vol = std::sqrt(total_var / T);
            }
        }

    public:
        // Generate fitted surface from the params
        static std::vector<VolSlice> GenerateFittedSurface(const std::vector<VolSlice>& market_surface, std::map<int, double>& forward_curve, FitterType type) {
            std::vector<VolSlice> fitted_surface = market_surface;

            for (auto& slice : fitted_surface) {
                double T = slice.expiry;

                for (auto& pt : slice.points) {
                    UpdatePointWithModel(pt, T, slice, type);
                }
            }
            std::sort(fitted_surface.begin(), fitted_surface.end(), [](const VolSlice& a, const VolSlice& b) {
                return a.expiry < b.expiry; // Expiries in ascending order
                });
            return fitted_surface;
        }

        static void ExportToCSV(const std::string& filepath, const std::vector<VolSlice>& fitted_surface) {
            
            std::ofstream file(filepath);
            
            if (!file.is_open()) {
                std::cerr << " [Evaluator Error] Unable to create file " << filepath << std::endl;
                return;
            }

            file << std::fixed << std::setprecision(6);

            std::set<double> all_strikes;
            for (const auto& slice : fitted_surface) {
                for (const auto& pt : slice.points) {
                    all_strikes.insert(pt.strike);
                }
            }

            // 2. Write Header : "Strike" then all expiries 
            file << "Strike";
            for (const auto& slice : fitted_surface) {
                file << "," << slice.expiry;
            }
            file << "\n";

            // 3. Build a map look for vol by [expiry][strike]
            // Structure : map< expiry, map< strike, implied_vol > >
            std::map<double, std::map<double, double>> vol_matrix;
            for (const auto& slice : fitted_surface) {
                for (const auto& pt : slice.points) {
                    vol_matrix[slice.expiry][pt.strike] = pt.implied_vol;
                }
            }

            // 4. Write each row
            for (double strike : all_strikes) {
				file << strike; // 1st column is the strike

                for (const auto& slice : fitted_surface) {
                    double expiry = slice.expiry;

                    if (vol_matrix[expiry].find(strike) != vol_matrix[expiry].end()) {
                        file << "," << vol_matrix[expiry][strike];
                    }
                    else {
                        file << ",";
                    }
                }
                file << "\n";
            }

            file.close();
            std::cout << " -> Success: Matrix of fitted vols exported to " << filepath << std::endl;
        }
};