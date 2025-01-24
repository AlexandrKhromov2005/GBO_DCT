#ifndef DCT_HPP
#define DCT_HPP

#include <array>
#include <algorithm>
#include <vector>

constexpr int DCT_SIZE = 8;

using Matrix8x8uc = std::array<std::array<unsigned char, DCT_SIZE>, DCT_SIZE>;
using Matrix8x8d = std::array<std::array<double, DCT_SIZE>, DCT_SIZE>;

struct DCTBlocks {
    std::vector<Matrix8x8d> blocks;
    int block_count_x;
    int block_count_y;
};

void dct_func(const Matrix8x8uc& block, Matrix8x8d& dct_block);
void rev_dct_func(Matrix8x8uc& block, const Matrix8x8d& dct_block);
void apply_dct_to_blocks(const std::vector<Matrix8x8uc>& input_blocks, DCTBlocks& dct_blocks);
void apply_rev_dct_to_blocks(const DCTBlocks& dct_blocks, std::vector<Matrix8x8uc>& output_blocks);
#endif // DCT_HPP