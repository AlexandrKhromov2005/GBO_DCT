#include <vector>
#include <array>
#include <random> // Добавлен заголовок для работы с генераторами случайных чисел
#include "population.hpp"

using VEC_POP = std::array<double, 22>;

class GBO {
public:
    int n; // Размер популяции
    int vec_size = 22; // Размер вектора
    int m; // Количество итераций
    int th; // Нижняя граница
    double pr; // Вероятность срабатывания leo

    const Matrix8x8d original_dct;
    const Matrix8x8uc original_block;

    unsigned char mode;

    double alpha, betta;

    double const epsilon = 0.05;

    std::vector<VEC_POP> population; // Популяция

    VEC_POP x; // Текущий вектор популяции
    VEC_POP x_next; // Вектор новой популяции
    VEC_POP x_best; // Лучший вектор
    VEC_POP x_worst; // Худший вектор
    VEC_POP x1, x2, x3;
    VEC_POP x_r1, x_r2, x_r3, x_r4, x_p; // Случайные векторы

    std::array<double, 22> f_values; 
    int best_ind, worst_ind;
    double best_f, worst_f;

    // Конструктор и деструктор
    GBO(int inp_n, int inp_m, int inp_th, double inp_pr, Matrix8x8d inp_original_dct,
    Matrix8x8uc inp_original_blocks, unsigned char inp_mode, int inp_best_ind, int inp_worst_ind, std::array<double, 22> inp_f_values, std::vector<VEC_POP> inp_population);
    ~GBO();

    // Методы
    void gbo(); // Функция оптимизации, совмещающая GSR и LEO
    double gsr(int cur_x); // Функция оптимизации GSR
    void leo(); // Функция "встряхивания" LEO

    double rand_0_1(); // Случайное число от 0 до 1
    double randn(); // Случайное число, распределенное по нормальному закону
    double rand_1_n(); // Случайное число от 1 до N
    double rand_m1_1(); // Случайное число от -1 до 1
    double dm(double& x_b_r1, double& x_n_r2); // Функция для вычисления dm

    void generate_population(); // Генерация популяции
    void rand_vecs(); // Генерация случайных векторов
    double new_rho(); // Вычисление нового значения rho

    void upgr_cur_x(int i);

private:
    // Генераторы случайных чисел
    std::random_device rd;
    std::mt19937 gen;
    std::uniform_real_distribution<> dis_0_1;
    std::normal_distribution<> dis_norm;
    std::uniform_int_distribution<> dis_1_n;
    std::uniform_real_distribution<> dis_m1_1;
};