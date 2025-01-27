#include <cmath>
#include <ctime>
#include <cstdlib>
#include <random>
#include <algorithm>
#include "gbo.hpp"

#include "population.hpp"

#define PI acos(-1)

// Конструктор
GBO::GBO(int inp_n, int inp_m, int inp_th, double inp_pr, Matrix8x8d inp_original_dct,
    Matrix8x8uc inp_original_blocks, unsigned char inp_mode)
    : mode(inp_mode), original_dct(inp_original_dct), original_block(inp_original_blocks), n(inp_n), m(inp_m), th(inp_th), pr(inp_pr),
      gen(rd()), dis_0_1(0.0, 1.0), dis_norm(0.0, 1.0),
      dis_1_n(1, n), dis_m1_1(-1.0, 1.0) {}

// Деструктор
GBO::~GBO() {}

// Генерация случайного числа от 0 до 1
double GBO::rand_0_1() {
    return dis_0_1(gen);
}

// Генерация случайного числа, распределенного по нормальному закону
double GBO::randn() {
    return dis_norm(gen);
}

// Генерация случайного числа от 1 до n
double GBO::rand_1_n() {
    return dis_1_n(gen);
}

// Генерация случайного числа от -1 до 1
double GBO::rand_m1_1() {
    return dis_m1_1(gen);
}

// Вычисление нового значения rho
double GBO::new_rho() {
    return (2 * rand_0_1() * alpha) - alpha;
}

// Генерация популяции
void GBO::generate_population() {
    for (int i = 0; i < n; i++) {
        VEC_POP temp_x;
        for (int j = 0; j < vec_size; j++) {
            temp_x[j] = -th + rand_0_1() * (th - (-th));
        }
        this->population.push_back(temp_x);
    }
}

// Функция оптимизации GSR
double GBO::gsr(int cur_x) {
    double delta = 2.0 * rand_0_1() * ((fabs(x_r1[cur_x] * x_r2[cur_x] * x_r3[cur_x] * x_r4[cur_x]) / 4.0) - x[cur_x]);
    double step = ((x_best[cur_x] - x_r1[cur_x]) + delta);
    double delta_x = rand_1_n() * fabs(step);
    double z = x[cur_x] - randn() * ((2 * delta_x * x[cur_x]) / (x_worst[cur_x] - x_best[cur_x] + epsilon));
    double p = rand_0_1() * (((z + x[cur_x]) / 2.0) + (rand_0_1() * delta_x));
    double q = rand_0_1() * (((z + x[cur_x]) / 2.0) - (rand_0_1() * delta_x));
    return randn() * new_rho() * ((2 * delta_x * x[cur_x]) / (p - q + epsilon));
}

// Функция для вычисления dm
double GBO::dm(double& x_b_r1, double& x_n_r2) {
    return rand_0_1() * new_rho() * (x_b_r1 - x_n_r2);
}

// Функция "встряхивания" LEO
void GBO::leo() {
    x = x_next;

    VEC_POP y = (rand_0_1() < 0.5) ? x : x_best;

    double mu1 = rand_0_1();
    double u1 = (mu1 < 0.5) ? 2 * rand_0_1() : 1;
    double u2 = (mu1 < 0.5) ? rand_0_1() : 1;
    double u3 = (mu1 < 0.5) ? rand_0_1() : 1;

    VEC_POP x_k;
    double mu2 = rand_0_1();
    if (mu2 < 0.5) {
        for (int j = 0; j < vec_size; j++) {
            x_k[j] = -th + rand_0_1() * (th - (-th));
        }
    } else {
        x_k = x_p;
    }

    for (int i = 0; i < vec_size; i++) {
        double f1 = rand_m1_1();
        double f2 = rand_m1_1();
        x_next[i] = y[i] + (f1 * (u1 * x_best[i] - u2 * x_k[i])) + (f2 * new_rho() * ((u3 * (x2[i] - x1[i]) + u2 * (x_r1[i] - x_r2[i])) / 2));
    }
}

// Основная функция оптимизации
void GBO::gbo() {
    XInd f;

    x_best = population[best_ind];
    best_f = f_values[best_ind];

    x_worst = population[worst_ind];
    worst_f = f_values[worst_ind];

    for (int cur_iter = 0; cur_iter < m; cur_iter++) {
        double betta = 0.2 + (1.2 - 0.2) * pow(1.0 - pow(((double)cur_iter / (double)m), 3), 2);
        alpha = fabs(betta * sin(((3.0 * PI) / 2.0) + sin((3 * PI * betta) / 2.0)));

        #pragma omp parallel for
        for (int cur_vec = 0; cur_vec < n; cur_vec++) {
            x = population[cur_vec];

            rand_vecs();

            for (int i = 0; i < vec_size; i++) {
                x1[i] = x[i] - gsr(i) + dm(x_best[i], x[i]);
                x2[i] = x_best[i] - gsr(i) + dm(x_r1[i], x_r2[i]);
                x3[i] = x[i] - new_rho() * (x2[i] - x1[i]);

                double r_a = rand_0_1();
                double r_b = rand_0_1();
                x_next[i] = r_a * (r_b * x1[i] + (1 - r_b) * x2[i]) + (1 - r_a) * x3[i];
            }

            if (rand_0_1() < pr) {
                leo();
            }
        
            upgr_cur_x(cur_vec);
        }
    }
}

// Генерация случайных векторов
void GBO::rand_vecs() {
    std::vector<int> ind_vec;

    while (ind_vec.size() < 4) {
        int ind = rand_1_n() - 1;
        auto it = std::find(ind_vec.begin(), ind_vec.end(), ind);

        if (it != ind_vec.end()) {
            continue;
        } else {
            ind_vec.push_back(ind);
        }
    }

    x_r1 = population[ind_vec[0]];
    x_r2 = population[ind_vec[1]];
    x_r3 = population[ind_vec[2]];
    x_r4 = population[ind_vec[3]];

    x_p = population[rand_1_n() - 1];
}

void GBO::upgr_cur_x(int i){
    PopulationOptimizer opt(th, n);

    Matrix8x8d modified_dct;
    // Применение преобразования X к DCT-коэффициентам
    opt.apply_x_transform(original_dct, x, modified_dct);
        
    // Обратное DCT-преобразование для получения пикселей
    Matrix8x8uc reconstructed_block;
    rev_dct_func(reconstructed_block, modified_dct);
        
    // Расчет фитнес-функции
    double fitness = opt.calculate_fitness(
        original_dct,
        modified_dct,
        original_block,
        mode
    );


    if(fitness < f_values[i]){
        population[i] = x;
        f_values[i] = fitness;
        if(fitness < best_f){
            x_best = x;
            best_f = fitness;
        }
    }
    else{
        x_worst = x;
        worst_f = fitness; 
    }
}