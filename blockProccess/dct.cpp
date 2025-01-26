#include "dct.hpp"

void DCTProcessor::dct_transform(const Block8x8& input, DoubleBlock8x8& output) {
    DoubleBlock8x8 double_block;
    convert_to_double(input, double_block);
    
    DoubleBlock8x8 temp;
    multiply_matrices(double_block, c_t, temp);
    multiply_matrices(c, temp, output);
}

void DCTProcessor::inverse_dct(const DoubleBlock8x8& input, Block8x8& output) {
    DoubleBlock8x8 temp;
    multiply_matrices(input, c, temp);
    
    DoubleBlock8x8 double_result;
    multiply_matrices(c_t, temp, double_result);
    convert_to_uchar(double_result, output);
}

void DCTProcessor::multiply_matrices(const DoubleBlock8x8& a, const DoubleBlock8x8& b, DoubleBlock8x8& result) {
    for(int i = 0; i < 8; ++i) {
        for(int j = 0; j < 8; ++j) {
            result[i][j] = 0.0;
            for(int k = 0; k < 8; ++k) {
                result[i][j] += a[i][k] * b[k][j];
            }
        }
    }
}

// Конвертация в double
void DCTProcessor::convert_to_double(const Block8x8& src, DoubleBlock8x8& dst) {
    for(int i = 0; i < 8; ++i) {
        for(int j = 0; j < 8; ++j) {
            dst[i][j] = static_cast<double>(src[i][j]);
        }
    }
}

// Конвертация в unsigned char с округлением и ограничением
void DCTProcessor::convert_to_uchar(const DoubleBlock8x8& src, Block8x8& dst) {
    for(int i = 0; i < 8; ++i) {
        for(int j = 0; j < 8; ++j) {
            double val = src[i][j];
            // Ручное округление
            val = (val >= 0.0) ? (val + 0.5) : (val - 0.5);
            // Ограничение диапазона и приведение типа
            dst[i][j] = static_cast<unsigned char>(
                std::clamp(static_cast<int>(val), 0, 255)
            );
        }
    }
}