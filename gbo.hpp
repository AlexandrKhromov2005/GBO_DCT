#include <vector>
#include <array>

using VEC_POP = std::array<double, 22>;

class GBO{
public:
    int n; // Размер популяции
    int m; // Количество итераций
    int th; // Нижняя граница
    double pr; //Вероятность срабатывания leo
    double cur_iter; // Текущая популяция

    double alpha, betta;

    double const epsilon = 0.05;

    std::vector<VEC_POP> population; // Популяция

    VEC_POP x; // Текущий вектор вектор
    VEC_POP x_best; // Лучший вектор
    VEC_POP x_worst; // Худший вектор
    VEC_POP x_r1, x_r2, x_r3, x_r4; // Случайные векторы


    void gbo(); // Функция оптимизации совмещающая GSR и LEO
    void gsr(int cur_x); // Функция оптимизации GSR 
    void leo(); // Функция "встряхивания" LEO

    double rand_0_1(); // Случайное число от 0 до 1 
    double randn(); // Случайных число, распределенное по нормальному закону
    double rand_1_n(); // Случайное число от 1 до N 

    void generate_population(); // Генерация популяции
    double new_rho();
};