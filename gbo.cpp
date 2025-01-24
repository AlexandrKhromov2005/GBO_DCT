#include <cmath>
#include <ctime>
#include <cstdlib>
#include <random> 

#include "gbo.hpp"

double GBO::rand_0_1(){
    return static_cast<double>(std::rand()) / RAND_MAX;
}

double randn() {
    // Инициализация генератора случайных чисел
    std::random_device rd;  // Источник энтропии
    std::mt19937 gen(rd()); // Генератор Mersenne Twister

    // Нормальное распределение с mean = 0 и stddev = 1
    std::normal_distribution<> d(0.0, 1.0);

    return d(gen);
}

double GBO::rand_1_n() {
    // Генерация случайного числа в диапазоне [0, RAND_MAX]
    double random_value = static_cast<double>(std::rand()) / RAND_MAX;
    // Масштабирование и сдвиг для получения диапазона [1, n]
    return random_value * n + 1;
}

double GBO::new_rho() {
	return (2 * rand_0_1() * alpha) - alpha;
}

void GBO::generate_population(){
    for(int i = 0; i < n; i++){
        VEC_POP x;
        for(int j = 0; j < 22; j++){
            x[j] = -th + rand_0_1() * (th - (-th));
        }
        this->population.push_back(x);
    }
}

void GBO::gsr(int cur_x){
    double delta = 2.0 * rand_0_1() * ((fabs(x_r1[cur_x] * x_r2[cur_x] * x_r3[cur_x] * x_r4[cur_x]) / 4.0) - x[cur_x]);
   
    double step = ((x_best[cur_x] - x_r1[cur_x]) + delta);

    double delta_x = rand_1_n() * fabs(step);

    double z = x[cur_iter] - randn() * ((2 * delta_x * x[cur_x]) / (x_worst[cur_x] - x_best[cur_x] + epsilon));

    double p = rand_0_1() * (((z + x[cur_x]) / 2.0) + (rand_0_1() * delta_x));
        
    double q = rand_0_1() * (((z + x[cur_x]) / 2.0) - (rand_0_1() * delta_x));

    double gsr = randn() * new_rho() * ((2 * delta_x * x[cur_x]) / (p - q + epsilon));

}