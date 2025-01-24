#include "population.hpp"
#include <algorithm>

PopulationOptimizer::PopulationOptimizer(double threshold, size_t population_size)
    : threshold_(threshold), population_size_(population_size),
      gen_(std::random_device{}()), dist_(-threshold, threshold) {
    initialize_population();
}

void PopulationOptimizer::initialize_population() {
    population_.resize(population_size_);
    for (auto& individual : population_) {
        for (auto& gene : individual) {
            gene = dist_(gen_);
        }
    }
}

XInd PopulationOptimizer::evaluate_population(
    const DCTBlocks& original_dct,
    const std::vector<Matrix8x8uc>& original_blocks,
    char mode
) {
    XInd results;
    double min_fitness = std::numeric_limits<double>::max();
    double max_fitness = -std::numeric_limits<double>::max();

    #pragma omp parallel for
    for (size_t i = 0; i < population_.size(); ++i) {
        Matrix8x8d modified_dct;
        apply_x_transform(original_dct.blocks[0], population_[i], modified_dct);
        
        Matrix8x8uc reconstructed_block;
        rev_dct_func(reconstructed_block, modified_dct);
        
        double fitness = calculate_fitness(
            original_dct.blocks[0],
            modified_dct,
            original_blocks[0],
            mode
        );

        #pragma omp critical
        {
            results.f_values[i] = fitness;
            if (fitness < min_fitness) {
                min_fitness = fitness;
                results.best = i;
            }
            if (fitness > max_fitness) {
                max_fitness = fitness;
                results.worst = i;
            }
        }
    }
    
    return results;
}

double PopulationOptimizer::calculate_fitness(
    const Matrix8x8d& original_dct,
    const Matrix8x8d& modified_dct,
    const Matrix8x8uc& original_block,
    char mode
) {
    Matrix8x8uc reconstructed_block;
    rev_dct_func(reconstructed_block, modified_dct);
    double mse = calculate_mse_block(original_block, reconstructed_block);
    double psnr = calculate_psnr_block(mse);

    double s0 = 0.0, s1 = 0.0;
    for (const auto& row : modified_dct) {
        for (double val : row) {
            s0 += std::abs(val);
            s1 += val * val;
        }
    }

    const double ratio = (mode == 0) ? (s1 / s0) : (s0 / s1);
    return ratio - 0.01 * psnr;
}

void PopulationOptimizer::apply_x_transform(
    const Matrix8x8d& src,
    const std::array<double, 22>& x,
    Matrix8x8d& dst
) {
    const std::array<std::pair<int, int>, 22> indices = {{
        {6,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5}, {0,6}, {0,7},
        {1,6}, {2,5}, {3,4}, {4,3}, {5,2}, {6,1}, {7,0}, {7,1},
        {6,2}, {5,3}, {4,4}, {3,5}, {2,6}, {1,7}
    }};

    dst = src;
    for (size_t i = 0; i < x.size(); ++i) {
        const auto& [row, col] = indices[i];
        dst[row][col] = sign(dst[row][col]) * std::abs(dst[row][col] + x[i]);
    }
}