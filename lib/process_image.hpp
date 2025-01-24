#pragma once
#include "dct.hpp"

struct ImageBlocks {
    std::vector<Matrix8x8uc> blocks;
    int block_count_x;
    int block_count_y;
};

unsigned char* crop_and_grayscale(
    const unsigned char* input,
    int width,
    int height,
    int channels,
    int* new_width,
    int* new_height
);

ImageBlocks split_into_blocks(const unsigned char* image, int width, int height);
unsigned char* assemble_image(const ImageBlocks& blocks, int* out_width, int* out_height);