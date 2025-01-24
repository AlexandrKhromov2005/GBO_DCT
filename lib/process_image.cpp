#include "process_image.hpp"
#include <cmath>

unsigned char* crop_and_grayscale(
    const unsigned char* input,
    int width,
    int height,
    int channels,
    int* new_width,
    int* new_height
) {
    *new_width = width - (width % 8);
    *new_height = height - (height % 8);
    if (*new_width == 0) *new_width = 8;
    if (*new_height == 0) *new_height = 8;

    unsigned char* output = new unsigned char[*new_width * *new_height];

    for (int y = 0; y < *new_height; ++y) {
        for (int x = 0; x < *new_width; ++x) {
            int input_idx = (y * width + x) * channels;
            unsigned char r = input[input_idx];
            unsigned char g = (channels > 1) ? input[input_idx + 1] : r;
            unsigned char b = (channels > 2) ? input[input_idx + 2] : r;
            output[y * *new_width + x] = static_cast<unsigned char>(0.299 * r + 0.587 * g + 0.114 * b);
        }
    }

    return output;
}

ImageBlocks split_into_blocks(const unsigned char* image, int width, int height) {
    ImageBlocks result;
    result.block_count_x = width / 8;
    result.block_count_y = height / 8;

    for (int by = 0; by < result.block_count_y; ++by) {
        for (int bx = 0; bx < result.block_count_x; ++bx) {
            Matrix8x8uc block;

            for (int y = 0; y < 8; ++y) {
                for (int x = 0; x < 8; ++x) {
                    int src_x = bx * 8 + x;
                    int src_y = by * 8 + y;
                    block[y][x] = image[src_y * width + src_x];
                }
            }

            result.blocks.push_back(block);
        }
    }

    return result;
}

unsigned char* assemble_image(const ImageBlocks& blocks, int* out_width, int* out_height) {
    *out_width = blocks.block_count_x * 8;
    *out_height = blocks.block_count_y * 8;

    unsigned char* image = new unsigned char[*out_width * *out_height];

    for (size_t i = 0; i < blocks.blocks.size(); ++i) {
        int bx = i % blocks.block_count_x;
        int by = i / blocks.block_count_x;
        const auto& block = blocks.blocks[i];

        for (int y = 0; y < 8; ++y) {
            for (int x = 0; x < 8; ++x) {
                int dst_x = bx * 8 + x;
                int dst_y = by * 8 + y;
                image[dst_y * *out_width + dst_x] = block[y][x];
            }
        }
    }

    return image;
}