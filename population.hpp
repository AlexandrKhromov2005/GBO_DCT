#ifndef POPULATION_HPP
#define POPULATION_HPP

#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <limits>
#include "dct.hpp"
#include "block_metrics.hpp"

struct XInd {
    int best = -1;
    int worst = -1;
    std::array<double, 22> f_values{};
};

class PopulationOptimizer {
public:
    PopulationOptimizer(double threshold, size_t population_size);

    void initialize_population();
    XInd evaluate_population(const DCTBlocks& original_dct, const std::vector<Matrix8x8uc>& original_blocks, char mode);
    
    // Геттеры
    const std::vector<std::array<double, 22>>& get_population() const { return population_; }
    
private:
    double threshold_;
    size_t population_size_;
    std::vector<std::array<double, 22>> population_;
    std::mt19937 gen_;
    std::uniform_real_distribution<double> dist_;

    double calculate_fitness(
        const Matrix8x8d& original_dct,
        const Matrix8x8d& modified_dct,
        const Matrix8x8uc& original_block,
        char mode
    );

    static void apply_x_transform(const Matrix8x8d& src, const std::array<double, 22>& x, Matrix8x8d& dst);
    static double sign(double val) { return (val >= 0) ? 1.0 : -1.0; }
};

#endif // POPULATION_HPP