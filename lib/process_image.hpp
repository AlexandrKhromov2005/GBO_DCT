#pragma once
#include "../blockProccess/dct.hpp"

#include <string>

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

ImageBlocks split_into_blocks(const std::vector<unsigned char>  image, int width, int height);
std::vector<unsigned char> assemble_image(const ImageBlocks& blocks, int* out_width, int* out_height);
std::vector<unsigned char> import_image(const std::string& filepath, int width, int height, int channels);
void export_image(const std::string& filepath, std::vector<unsigned char> image_vec, int width, int height, int channels);