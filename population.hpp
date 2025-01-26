#pragma once
#include <array>
#include <vector>
#include <limits>
#include <cmath>
#include <cstdlib>
#include "block_metrics.hpp"
#include "dct.hpp"

// Структура для хранения индексов лучшего и худшего решения
struct XIndices {
    int best = -1;          // Индекс лучшей особи в популяции
    int worst = -1;         // Индекс худшей особи в популяции
    std::vector<double> f_values; // Значения целевой функции для всех особей
};

// Класс для работы с популяцией решений
class PopulationOptimizer {
public:
    using Block8x8 = DCTProcessor::Block8x8;
    using DoubleBlock8x8 = DCTProcessor::DoubleBlock8x8;
    
    // Конструктор с параметрами
    PopulationOptimizer(double threshold, int population_size);
    
    // Создание начальной популяции
    void create_population(std::vector<std::array<double, 22>>& population);
    
    // Поиск лучшего решения в популяции
    int find_best(const std::vector<std::array<double, 22>>& population,
                 const DoubleBlock8x8& original_dct,
                 const Block8x8& original_block,
                 char bit) const;
    
    // Поиск лучшего и худшего решений
    XIndices find_bw(const std::vector<std::array<double, 22>>& population,
                   const DoubleBlock8x8& original_dct,
                   const Block8x8& original_block,
                   char bit) const;

private:
    double threshold_;       // Порог для генерации начальных значений
    int population_size_;    // Размер популяции
    
    // Генерация случайного числа [0, 1)
    static double rand_double();
    
    // Целевая функция для оценки решения
    double objective_function(const DoubleBlock8x8& original_dct,
                             const Block8x8& original_block,
                             const std::array<double, 22>& x,
                             char bit) const;
    
    // Применение модификаций к блоку DCT
    void apply_x(const DoubleBlock8x8& original,
                const std::array<double, 22>& x,
                DoubleBlock8x8& modified) const;
    
    // Вспомогательные функции для вычисления S0 и S1
    static double calculate_s0(const DoubleBlock8x8& block);
    static double calculate_s1(const DoubleBlock8x8& block);
};