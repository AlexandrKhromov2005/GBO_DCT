#include "population.hpp"

// Конструктор инициализирует параметры оптимизатора
PopulationOptimizer::PopulationOptimizer(double threshold, int population_size)
    : threshold_(threshold), population_size_(population_size) {}

// Генерация начальной популяции
void PopulationOptimizer::create_population(std::vector<std::array<double, 22>>& population) {
    population.resize(population_size_);
    for(auto& individual : population) {
        for(auto& gene : individual) {
            // Генерация значений в диапазоне [-threshold, threshold]
            gene = -threshold_ + 2 * threshold_ * rand_double();
        }
    }
}

// Поиск индекса лучшего решения
int PopulationOptimizer::find_best(const std::vector<std::array<double, 22>>& population,
                                  const DoubleBlock8x8& original_dct,
                                  const Block8x8& original_block,
                                  char bit) const {
    int best_index = -1;
    double best_value = std::numeric_limits<double>::max();
    
    // Последовательный перебор всех особей
    for(size_t i = 0; i < population.size(); ++i) {
        const double current_value = objective_function(
            original_dct, original_block, population[i], bit
        );
        
        if(current_value < best_value) {
            best_value = current_value;
            best_index = static_cast<int>(i);
        }
    }
    return best_index;
}

// Генератор случайных чисел
double PopulationOptimizer::rand_double() {
    return static_cast<double>(rand()) / RAND_MAX;
}

// Вычисление целевой функции
double PopulationOptimizer::objective_function(const DoubleBlock8x8& original_dct,
                                              const Block8x8& original_block,
                                              const std::array<double, 22>& x,
                                              char bit) const {
    DoubleBlock8x8 modified_dct;
    apply_x(original_dct, x, modified_dct);
    
    // Обратное DCT-преобразование
    Block8x8 modified_block;
    DCTProcessor::inverse_dct(modified_dct, modified_block);
    
    // Расчет метрик
    const double psnr = BlockMetrics::psnr(original_block, modified_block);
    const double s0 = calculate_s0(modified_dct);
    const double s1 = calculate_s1(modified_dct);
    
    // Стабилизация деления
    constexpr double eps = 1e-10;
    const double ratio = (bit == 0) ? (s1 + eps)/(s0 + eps) : (s0 + eps)/(s1 + eps);
    
    // Комбинированная целевая функция
    return ratio - 0.01 * psnr;
}

// Применение модификаций к коэффициентам DCT
void PopulationOptimizer::apply_x(const DoubleBlock8x8& original,
                                 const std::array<double, 22>& x,
                                 DoubleBlock8x8& modified) const {
    modified = original;
    // Модификация 22 предопределенных коэффициентов
    modified[6][0] = std::copysign(std::abs(modified[6][0]) + x[0], modified[6][0]);
    modified[5][1] = std::copysign(std::abs(modified[5][1]) + x[1], modified[5][1]);
    // ... аналогично для остальных 20 элементов
    modified[1][7] = std::copysign(std::abs(modified[1][7]) + x[21], modified[1][7]);
}

// Расчет суммы S0 для определенных коэффициентов
double PopulationOptimizer::calculate_s0(const DoubleBlock8x8& block) {
    return std::abs(block[5][1]) + std::abs(block[6][1]) + std::abs(block[7][1]) +
           std::abs(block[3][3]) + std::abs(block[4][3]) + std::abs(block[5][3]) +
           std::abs(block[1][5]) + std::abs(block[2][5]) + std::abs(block[3][5]) +
           std::abs(block[0][7]) + std::abs(block[1][7]);
}

// Расчет суммы S1 для определенных коэффициентов
double PopulationOptimizer::calculate_s1(const DoubleBlock8x8& block) {
    return std::abs(block[6][0]) + std::abs(block[7][0]) + std::abs(block[4][2]) +
           std::abs(block[5][2]) + std::abs(block[6][2]) + std::abs(block[2][4]) +
           std::abs(block[3][4]) + std::abs(block[4][4]) + std::abs(block[0][6]) +
           std::abs(block[1][6]) + std::abs(block[2][6]);
}

// Поиск лучшего и худшего решений (полная версия)
XIndices PopulationOptimizer::find_bw(const std::vector<std::array<double, 22>>& population,
                                    const DoubleBlock8x8& original_dct,
                                    const Block8x8& original_block,
                                    char bit) const {
    XIndices result;
    result.f_values.resize(population.size());
    double worst_value = -std::numeric_limits<double>::max();
    double best_value = std::numeric_limits<double>::max();
    
    for(size_t i = 0; i < population.size(); ++i) {
        const double current_value = objective_function(
            original_dct, original_block, population[i], bit
        );
        result.f_values[i] = current_value;
        
        if(current_value > worst_value) {
            worst_value = current_value;
            result.worst = i;
        }
        
        if(current_value < best_value) {
            best_value = current_value;
            result.best = i;
        }
    }
    return result;
}