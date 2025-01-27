#include <ctime>
#include <cstdlib>
#include <iostream>
#include "gbo.hpp"
#include "population.hpp"
#include "lib/process_image.hpp"
#include "blockProccess/dct.hpp"
#include "blockProccess/block_metrics.hpp"
#include <chrono> // Для измерения времени


int main(){
    // Засекаем время начала выполнения
    auto start = std::chrono::high_resolution_clock::now();

    // Инициализация генератора случайных чисел
    std::srand(std::time(0)); 

    std::vector<unsigned char> img_pixels = import_image("images/lenna.png", 512, 512, 1);

    ImageBlocks image = split_into_blocks(img_pixels, 512, 512);

    for(int i = 0; i < image.blocks.size(); i++){
        std::cout << i << "\n"; 

        PopulationOptimizer population(10, 30);
        
        Matrix8x8d dct_block;
        dct_func(image.blocks[i], dct_block); 

        XInd f_result = population.evaluate_population(dct_block, image.blocks[i], 0);

        GBO optimizer(30, 40, 10, 0.5, dct_block, image.blocks[i], 0);

        optimizer.population = population.get_population();

        optimizer.best_ind = f_result.best;
        optimizer.worst_ind = f_result.worst;
        optimizer.f_values = f_result.f_values;



        optimizer.gbo();
    }

    PopulationOptimizer population(10, 30);
    XInd pop_result;

    //pop_result = population.evaluate_population();


    // Засекаем время окончания выполнения
    auto end = std::chrono::high_resolution_clock::now();

    // Вычисляем разницу между началом и концом в секундах
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    // Выводим время выполнения в секундах
    std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;

    return 0; 
}