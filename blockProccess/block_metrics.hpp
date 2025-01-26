#pragma once
#include <array>

class BlockMetrics {
public:
    using Block8x8 = std::array<std::array<unsigned char, 8>, 8>;
    
    static double mse(const Block8x8& original, const Block8x8& changed);
    static double psnr(const Block8x8& original, const Block8x8& changed);
};