#include "population.hpp"
#include <algorithm>

// Конструктор
PopulationOptimizer::PopulationOptimizer(double threshold, size_t population_size)
    : threshold_(threshold), population_size_(population_size),
      gen_(std::random_device{}()), dist_(-threshold, threshold) {
    initialize_population();
}

// Инициализация популяции
void PopulationOptimizer::initialize_population() {
    population_.resize(population_size_);
    for (auto& individual : population_) {
        for (auto& gene : individual) {
            gene = dist_(gen_);
        }
    }
}

// Оценка популяции
XInd PopulationOptimizer::evaluate_population(
    const Matrix8x8d& original_dct,
    const Matrix8x8uc& original_block,
    char mode
) {
    XInd results;
    double min_fitness = std::numeric_limits<double>::max();
    double max_fitness = -std::numeric_limits<double>::max();

    #pragma omp parallel for
    for (size_t i = 0; i < population_.size(); ++i) {
        Matrix8x8d modified_dct;
        apply_x_transform(original_dct, population_[i], modified_dct);
        
        Matrix8x8uc reconstructed_block;
        rev_dct_func(reconstructed_block, modified_dct);
        
        double fitness = calculate_fitness(
            original_dct,
            modified_dct,
            original_block,
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

// Расчет фитнес-функции
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

    // Индексы для 1-region (S1) и 0-region (S0)
    const std::array<std::pair<int, int>, 11> s1_indices = {{
        {0,7}, {0, 6}, {2, 6}, {2, 5}, {2, 4}, {4, 4}, {4, 3}, {4, 2}, {6, 2}, {6, 1}, {6, 0}
    }};
    
    const std::array<std::pair<int, int>, 11> s0_indices = {{
        {1,7}, {1, 6}, {1, 5}, {3, 5}, {3, 4}, {3, 3}, {5, 3}, {5, 2}, {5, 1}, {7, 1}, {7, 0}
    }};

    // Расчет S1 и S0 для соответствующих регионов
    double s1 = 0.0, s0 = 0.0;
    
    for (const auto& [row, col] : s1_indices) {
        double val = modified_dct[row][col];
        s1 += std::abs(val);  // Сумма квадратов для 1-region
    }
    
    
    for (const auto& [row, col] : s0_indices) {
        double val = modified_dct[row][col];
        s0 += std::abs(val);  // Сумма модулей для 0-region
    }

    if (mode == 0) {
        if (s0 == 0){ s0 = 1-e9;}
            const double ratio = s1/s0;
    } else {
        if (s1 == 0){ s1 = 1-e9;}
            const double ratio = s0/s1;
    }
    return ratio - 0.01 * psnr;
}

// Применение X-преобразования
void PopulationOptimizer::apply_x_transform(
    const Matrix8x8d& src,
    const std::array<double, 22>& x,
    Matrix8x8d& dst
) {
    // Объединенные индексы (порядок важен для применения изменений)
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
