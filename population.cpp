#include "population.hpp"
#include <algorithm>

// Конструктор: инициализация генератора, распределения и популяции
PopulationOptimizer::PopulationOptimizer(double threshold, size_t population_size)
    : threshold_(threshold), population_size_(population_size),
      gen_(std::random_device{}()), dist_(-threshold, threshold) {
    initialize_population();
}

// Инициализация популяции случайными векторами изменений
void PopulationOptimizer::initialize_population() {
    population_.resize(population_size_);
    for (auto& individual : population_) {
        for (auto& gene : individual) {
            gene = dist_(gen_); // Генерация значений в диапазоне [-threshold, threshold]
        }
    }
}

// Оценка популяции: параллельный расчет фитнес-функции для каждого индивида
XInd PopulationOptimizer::evaluate_population(
    const Matrix8x8d& original_dct,
    const Matrix8x8uc& original_block,
    char mode
) {
    XInd results;
    double min_fitness = std::numeric_limits<double>::max();
    double max_fitness = -std::numeric_limits<double>::max();

    // Параллельный цикл для оценки каждой особи
    #pragma omp parallel for
    for (size_t i = 0; i < population_.size(); ++i) {
        Matrix8x8d modified_dct;
        // Применение преобразования X к DCT-коэффициентам
        apply_x_transform(original_dct, population_[i], modified_dct);
        
        // Обратное DCT-преобразование для получения пикселей
        Matrix8x8uc reconstructed_block;
        rev_dct_func(reconstructed_block, modified_dct);
        
        // Расчет фитнес-функции
        double fitness = calculate_fitness(
            original_dct,
            modified_dct,
            original_block,
            mode
        );

        // Обновление лучших/худших результатов (критическая секция для потокобезопасности)
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

// // Оценка популяции: параллельный расчет фитнес-функции для каждого индивида
// XInd PopulationOptimizer::evaluate_population(
//     const DCTBlocks& original_dct,
//     const std::vector<Matrix8x8uc>& original_blocks,
//     char mode
// ) {
//     XInd results;
//     double min_fitness = std::numeric_limits<double>::max();
//     double max_fitness = -std::numeric_limits<double>::max();

//     // Параллельный цикл для оценки каждой особи
//     #pragma omp parallel for
//     for (size_t i = 0; i < population_.size(); ++i) {
//         Matrix8x8d modified_dct;
//         // Применение преобразования X к DCT-коэффициентам
//         apply_x_transform(original_dct.blocks[0], population_[i], modified_dct);
        
//         // Обратное DCT-преобразование для получения пикселей
//         Matrix8x8uc reconstructed_block;
//         rev_dct_func(reconstructed_block, modified_dct);
        
//         // Расчет фитнес-функции
//         double fitness = calculate_fitness(
//             original_dct.blocks[0],
//             modified_dct,
//             original_blocks[0],
//             mode
//         );

//         // Обновление лучших/худших результатов (критическая секция для потокобезопасности)
//         #pragma omp critical
//         {
//             results.f_values[i] = fitness;
//             if (fitness < min_fitness) {
//                 min_fitness = fitness;
//                 results.best = i;
//             }
//             if (fitness > max_fitness) {
//                 max_fitness = fitness;
//                 results.worst = i;
//             }
//         }
//     }
    
//     return results;
// }

// Расчет фитнес-функции:
// F = (S1/S0 или S0/S1) - 0.01 * PSNR (см. формулу (24) из статьи)
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

    // Индексы среднечастотных коэффициентов (совпадают с apply_x_transform)
    const std::array<std::pair<int, int>, 22> indices = {{
        {6,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5}, {0,6}, {0,7},
        {1,6}, {2,5}, {3,4}, {4,3}, {5,2}, {6,1}, {7,0}, {7,1},
        {6,2}, {5,3}, {4,4}, {3,5}, {2,6}, {1,7}
    }};

    double s0 = 0.0, s1 = 0.0;
    for (const auto& [row, col] : indices) {
        double val = modified_dct[row][col];
        s0 += std::abs(val);  // Сумма модулей среднечастотных коэффициентов
        s1 += val * val;      // Сумма квадратов среднечастотных коэффициентов
    }

    const double ratio = (mode == 0) ? (s1 / s0) : (s0 / s1);
    return ratio - 0.01 * psnr;
}

// Применение преобразования X к коэффициентам DCT:
// Изменяет коэффициенты в позициях, указанных в indices (см. рис. 1b из статьи)
void PopulationOptimizer::apply_x_transform(
    const Matrix8x8d& src,
    const std::array<double, 22>& x,
    Matrix8x8d& dst
) {
    // Индексы коэффициентов в среднечастотной области (22 позиции)
    const std::array<std::pair<int, int>, 22> indices = {{
        {6,0}, {5,1}, {4,2}, {3,3}, {2,4}, {1,5}, {0,6}, {0,7},
        {1,6}, {2,5}, {3,4}, {4,3}, {5,2}, {6,1}, {7,0}, {7,1},
        {6,2}, {5,3}, {4,4}, {3,5}, {2,6}, {1,7}
    }};

    dst = src;
    for (size_t i = 0; i < x.size(); ++i) {
        const auto& [row, col] = indices[i];
        // Изменение коэффициента с сохранением знака (формула (23) из статьи)
        dst[row][col] = sign(dst[row][col]) * std::abs(dst[row][col] + x[i]);
    }
}
