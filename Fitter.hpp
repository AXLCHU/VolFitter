#pragma once
#include "Models.hpp"
#include <vector>

class VolFitter {
public:
    static void FitSurface(std::vector<VolSlice>& slices, FitterType type);

private:
    static SVIParameters FitSingleSliceSVI(const VolSlice& slice, const SVIParameters* prev_params = nullptr);
    static SSVIParameters FitSingleSliceSSVI(const VolSlice& slice, const SSVIParameters* prev_params, double penalty_scaler);
    static double DeduceAtmVol(const VolSlice& slice);
    static SVIJWParameters FitSingleSliceSVIJW(const VolSlice& slice, const SVIJWParameters* prev_params, double penalty_scaler);
};
