#define _USE_MATH_DEFINES 
#include "MarketDataReader.hpp"
#include <fstream>
#include <sstream>
#include <cmath>
#include <iostream>
#include <stdexcept>
#include <algorithm>

// Fonction de parsing des dates (Format JJ/MM/AAAA) renvoyant le nombre absolu de jours cumulés
long GetAbsoluteDays(int d, int m, int y) {
    if (m < 3) { y--; m += 12; }
    return 365 * y + y / 4 - y / 100 + y / 400 + (153 * m + 2) / 5 + d;
}

// Calcule l'écart en jours entre la date cible et la date d'évaluation
int GetDaysToMaturity(const std::string& date_str, int ref_D = 30, int ref_M = 6, int ref_Y = 2025) {
    std::stringstream ss(date_str);
    std::string item;
    int day = 1, month = 1, year = 2025;
    try {
        if (std::getline(ss, item, '/')) day = std::stoi(item);
        if (std::getline(ss, item, '/')) month = std::stoi(item);
        if (std::getline(ss, item, '/')) year = std::stoi(item);
    }
    catch (...) { return 0; }

    return static_cast<int>(GetAbsoluteDays(day, month, year) - GetAbsoluteDays(ref_D, ref_M, ref_Y));
}

double CalculateBlackScholesVega(double strike, double forward, double expiry, double vol) {
    if (expiry <= 0.0 || vol <= 0.0) return 0.01;
    double std_dev = vol * std::sqrt(expiry);
    double d1 = (std::log(forward / strike) + 0.5 * std_dev * std_dev) / std_dev;
    return forward * std::sqrt(expiry) * (std::exp(-0.5 * d1 * d1) / std::sqrt(2.0 * M_PI));
}

// 1. LECTURE DU FICHIER FORWARD
std::map<int, double> MarketDataReader::ReadForwardCurve(const std::string& filepath) {
    std::map<int, double> forward_curve;
    std::ifstream file(filepath);
    if (!file.is_open()) throw std::runtime_error("Impossible d'ouvrir le fichier forward : " + filepath);

    std::string line;
    std::getline(file, line); // Ignorer l'en-tête (days,forward,rate)

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream ss(line);
        std::string cell;

        std::getline(ss, cell, ';'); int days = std::stoi(cell);
        std::getline(ss, cell, ';'); double forward = std::stod(cell);

        forward_curve[days] = forward;
    }
    file.close();
    return forward_curve;
}

// 2. LECTURE ET ALIGNEMENT DES VOLATILITÉS AVEC LES FORWARDS DYNAMIQUES
std::vector<VolSlice> MarketDataReader::ReadCSV(const std::string& filepath, const std::map<int, double>& forward_curve) {
    std::vector<VolSlice> slices;
    std::ifstream file(filepath);
    if (!file.is_open()) throw std::runtime_error("Volatility file is locked : " + filepath);

    std::string line;
    std::getline(file, line); // Ligne d'en-tête des dates

    std::stringstream header_stream(line);
    std::string cell;
    std::getline(header_stream, cell, ';'); // Ignorer la cellule "Strike"

    // Configuration de la date d'évaluation (à ajuster selon la date d'extraction de vos fichiers)
    // Ici configuré en juin 2025 pour correspondre de manière cohérente aux échéances de juillet 2025 lues précédemment
    int ref_D = 30, ref_M = 6, ref_Y = 2025;

    while (std::getline(header_stream, cell, ';')) {
        while (!cell.empty() && (cell.back() == '\r' || cell.back() == '\n' || cell.back() == ' ')) cell.pop_back();
        if (!cell.empty()) {
            VolSlice slice;
            slice.expiry_days = GetDaysToMaturity(cell, ref_D, ref_M, ref_Y);
            slice.expiry = static_cast<double>(slice.expiry_days) / 365.25;
            slices.push_back(slice);
        }
    }

    while (std::getline(file, line)) {
        if (line.empty()) continue;
        std::stringstream line_stream(line);
        std::getline(line_stream, cell, ';');
        if (cell.empty()) continue;
        double strike = std::stod(cell);

        size_t col_idx = 0;
        while (std::getline(line_stream, cell, ';')) {
            if (col_idx >= slices.size()) break;
            if (!cell.empty() && cell != "\r" && cell != "\n") {
                double iv = std::stod(cell);
                double t = slices[col_idx].expiry;
                int days = slices[col_idx].expiry_days;

                // Récupération automatique du Forward correspondant aux "days" ou choix du plus proche
                double forward = 6000.0; // Valeur refuge par défaut
                if (forward_curve.count(days)) {
                    forward = forward_curve.at(days);
                }
                else {
                    // Si écart minime de jours, on cherche la clé la plus proche
                    auto it = forward_curve.lower_bound(days);
                    if (it != forward_curve.end()) forward = it->second;
                }

                double k = std::log(strike / forward); // k calculé avec le bon Forward de la slice !
                double w = iv * iv * t;
                double vega = CalculateBlackScholesVega(strike, forward, t, iv);

                slices[col_idx].points.push_back(MarketPoint{ strike, iv, vega, k, w });
            }
            col_idx++;
        }
    }
    file.close();
    return slices;
}
