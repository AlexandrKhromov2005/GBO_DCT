#include "dct.hpp"
#include <array>
#include <cmath>

namespace {
    constexpr Matrix8x8d c = {{
        {0.353553, 0.353553, 0.353553, 0.353553, 0.353553, 0.353553, 0.353553, 0.353553},
        {0.490393, 0.415735, 0.277785, 0.0975452, -0.0975452, -0.277785, -0.415735, -0.490393},
        {0.46194, 0.191342, -0.191342, -0.46194, -0.46194, -0.191342, 0.191342, 0.46194},
        {0.415735, -0.0975452, -0.490393, -0.277785, 0.277785, 0.490393, 0.0975452, -0.415735},
        {0.353553, -0.353553, -0.353553, 0.353553, 0.353553, -0.353553, -0.353553, 0.353553},
        {0.277785, -0.490393, 0.0975452, 0.415735, -0.415735, -0.0975452, 0.490393, -0.277785},
        {0.191342, -0.46194, 0.46194, -0.191342, -0.191342, 0.46194, -0.46194, 0.191342},
        {0.0975452, -0.277785, 0.415735, -0.490393, 0.490393, -0.415735, 0.277785, -0.0975452}
    }};

    constexpr Matrix8x8d c_t = {{
        {0.353553, 0.490393, 0.46194, 0.415735, 0.353553, 0.277785, 0.191342, 0.0975452},
        {0.353553, 0.415735, 0.191342, -0.0975452, -0.353553, -0.490393, -0.46194, -0.277785},
        {0.353553, 0.277785, -0.191342, -0.490393, -0.353553, 0.0975452, 0.46194, 0.415735},
        {0.353553, 0.0975452, -0.46194, -0.277785, 0.353553, 0.415735, -0.191342, -0.490393},
        {0.353553, -0.0975452, -0.46194, 0.277785, 0.353553, -0.415735, -0.191342, 0.490393},
        {0.353553, -0.277785, -0.191342, 0.490393, -0.353553, -0.0975452, 0.46194, -0.415735},
        {0.353553, -0.415735, 0.191342, 0.0975452, -0.353553, 0.490393, -0.46194, 0.277785},
        {0.353553, -0.490393, 0.46194, -0.415735, 0.353553, -0.277785, 0.191342, -0.0975452}
    }};

    inline void multiply_matrices(const Matrix8x8d& a, const Matrix8x8d& b, Matrix8x8d& result) noexcept {
        for (int i = 0; i < DCT_SIZE; ++i) {
            for (int j = 0; j < DCT_SIZE; ++j) {
                double sum = 0.0;
                for (int k = 0; k < DCT_SIZE; ++k) {
                    sum += a[i][k] * b[k][j];
                }
                result[i][j] = sum;
            }
        }
    }

    inline void convert_to_double(const Matrix8x8uc& src, Matrix8x8d& dest) noexcept {
        for (int i = 0; i < DCT_SIZE; ++i) {
            for (int j = 0; j < DCT_SIZE; ++j) {
                dest[i][j] = static_cast<double>(src[i][j]) - 128.0;
            }
        }
    }

    inline void convert_to_uchar(const Matrix8x8d& src, Matrix8x8uc& dest) noexcept {
        for (int i = 0; i < DCT_SIZE; ++i) {
            for (int j = 0; j < DCT_SIZE; ++j) {
                double val = src[i][j] + 128.0;
                val = std::clamp(val, 0.0, 255.0);
                dest[i][j] = static_cast<unsigned char>(val + 0.5);
            }
        }
    }
}

void dct_func(const Matrix8x8uc& block, Matrix8x8d& dct_block) {
    Matrix8x8d double_block;
    Matrix8x8d temp_block;
    convert_to_double(block, double_block);
    multiply_matrices(double_block, c_t, temp_block);
    multiply_matrices(c, temp_block, dct_block);
}

void rev_dct_func(Matrix8x8uc& block, const Matrix8x8d& dct_block) {
    Matrix8x8d temp_block;
    Matrix8x8d double_block;
    multiply_matrices(dct_block, c, temp_block);
    multiply_matrices(c_t, temp_block, double_block);
    convert_to_uchar(double_block, block);
}

void apply_dct_to_blocks(const std::vector<Matrix8x8uc>& input_blocks, DCTBlocks& dct_blocks) {
    dct_blocks.blocks.clear();
    dct_blocks.block_count_x = 0;
    dct_blocks.block_count_y = 0;
    dct_blocks.blocks.reserve(input_blocks.size());

    for (const auto& block : input_blocks) {
        Matrix8x8d dct_block;
        dct_func(block, dct_block);
        dct_blocks.blocks.push_back(dct_block);
    }
}

void apply_rev_dct_to_blocks(const DCTBlocks& dct_blocks, std::vector<Matrix8x8uc>& output_blocks) {
    output_blocks.clear();
    output_blocks.reserve(dct_blocks.blocks.size());

    for (const auto& dct_block : dct_blocks.blocks) {
        Matrix8x8uc output_block;
        rev_dct_func(output_block, dct_block);
        output_blocks.push_back(output_block);
    }
}
