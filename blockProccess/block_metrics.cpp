#include "block_metrics.hpp"

double calculate_mse_block(const Matrix8x8uc& original, const Matrix8x8uc& reconstructed) {
    double sum_squared_error = 0.0;
    
    for (int y = 0; y < DCT_SIZE; ++y) {
        for (int x = 0; x < DCT_SIZE; ++x) {
            double diff = static_cast<double>(original[y][x]) - 
                          static_cast<double>(reconstructed[y][x]);
            sum_squared_error += diff * diff;
        }
    }
    
    return sum_squared_error / (DCT_SIZE * DCT_SIZE); 
}

double calculate_psnr_block(double mse) {
    if (mse <= 0.0) {
        return INFINITY; 
    }
    const double max_pixel_value = 255.0;
    return 10.0 * std::log10((max_pixel_value * max_pixel_value) / mse);
}