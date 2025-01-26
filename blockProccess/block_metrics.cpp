#include "block_metrics.hpp"
#include <cmath>

double BlockMetrics::mse(const Block8x8& original, const Block8x8& changed) {
    double sum = 0.0;
    for(size_t i = 0; i < 8; ++i) {
        for(size_t j = 0; j < 8; ++j) {
            double diff = original[i][j] - changed[i][j];
            sum += diff * diff;
        }
    }
    return sum / 64.0;
}

double BlockMetrics::psnr(const Block8x8& original, const Block8x8& changed) {
    const double mse_val = mse(original, changed);
    return mse_val <= 1e-10 ? 100.0 : 10.0 * log10(65025.0 / mse_val);
}