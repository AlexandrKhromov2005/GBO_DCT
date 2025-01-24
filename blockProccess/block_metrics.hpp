#ifndef BLOCK_METRICS_HPP
#define BLOCK_METRICS_HPP

#include <cmath>
#include "dct.hpp"

double calculate_mse_block(const Matrix8x8uc& original, const Matrix8x8uc& reconstructed);
double calculate_psnr_block(double mse);

#endif BLOCK_METRICS_HPP
