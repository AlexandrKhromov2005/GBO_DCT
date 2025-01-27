#ifndef POPULATION_HPP
#define POPULATION_HPP

#include <array>
#include <vector>
#include <random>
#include <cmath>
#include <omp.h>
#include <limits>
#include "blockProccess/dct.hpp"
#include "blockProccess/block_metrics.hpp"

// Структура для хранения результатов оценки популяции:
// - best: индекс лучшей особи (с минимальным значением фитнес-функции)
// - worst: индекс худшей особи (с максимальным значением фитнес-функции)
// - f_values: значения фитнес-функции для каждой особи популяции
struct XInd {
    int best = -1;
    int worst = -1;
    std::array<double, 22> f_values{};
};

class PopulationOptimizer {
public:
    // Конструктор:
    // - threshold: порог изменения коэффициентов DCT (параметр Th из статьи)
    // - population_size: размер популяции для оптимизации (N)
    PopulationOptimizer(double threshold, size_t population_size);

    // Инициализация популяции случайными значениями в диапазоне [-threshold, threshold]
    void initialize_population();
    
    // Оценка популяции:
    // - original_dct: исходные DCT-блоки
    // - original_blocks: исходные пиксельные блоки
    // - mode: режим (0 или 1), определяющий соотношение сумм коэффициентов (S1/S0 или S0/S1)
    // Возвращает структуру XInd с результатами оценки
    XInd evaluate_population(const Matrix8x8d& original_dct, const Matrix8x8uc& original_block, char mode);
    
    // Геттер для получения текущей популяции
    const std::vector<std::array<double, 22>>& get_population() const { return population_; }

    // Сочетает PSNR (для незаметности) и соотношение сумм коэффициентов (для устойчивости)
    double calculate_fitness(
        const Matrix8x8d& original_dct,
        const Matrix8x8d& modified_dct,
        const Matrix8x8uc& original_block,
        char mode
    );

    // Применение преобразования X к DCT-коэффициентам:
    // Изменяет коэффициенты в соответствии с вектором x (см. статью, раздел 3.1)
    static void apply_x_transform(const Matrix8x8d& src, const std::array<double, 22>& x, Matrix8x8d& dst);
    
private:
    double threshold_;          // Порог изменения коэффициентов (Th)
    size_t population_size_;    // Размер популяции
    std::vector<std::array<double, 22>> population_; // Популяция векторов изменений (22 коэффициента)
    std::mt19937 gen_;          // Генератор случайных чисел
    std::uniform_real_distribution<double> dist_; // Распределение для инициализации

    // Расчет фитнес-функции:
    // // Сочетает PSNR (для незаметности) и соотношение сумм коэффициентов (для устойчивости)
    // double calculate_fitness(
    //     const Matrix8x8d& original_dct,
    //     const Matrix8x8d& modified_dct,
    //     const Matrix8x8uc& original_block,
    //     char mode
    // );

    // // Применение преобразования X к DCT-коэффициентам:
    // // Изменяет коэффициенты в соответствии с вектором x (см. статью, раздел 3.1)
    // static void apply_x_transform(const Matrix8x8d& src, const std::array<double, 22>& x, Matrix8x8d& dst);
    
    // Вспомогательная функция для определения знака
    static double sign(double val) { return (val >= 0) ? 1.0 : -1.0; }
};

#endif // POPULATION_HPP
