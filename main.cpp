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

    std::vector<unsigned char> wm =  import_image("images/wm.png", 32, 32, 1);
    for(int i = 0; i < wm.size(); i++){
        if(wm[i] > 127){ wm[i] = 1;}
        else{ wm[i] = 0;}
    }

    for(int i = 0; i < image.blocks.size(); i++){
        std::cout << i << "\n"; 

        PopulationOptimizer population(10, 30);
        
        Matrix8x8d dct_block;
        dct_func(image.blocks[i], dct_block); 

        XInd f_result = population.evaluate_population(dct_block, image.blocks[i], wm[i % wm.size()]);

        // GBO optimizer(30, 40, 10, 0.5, dct_block, image.blocks[i], wm[i % wm.size()], f_result.best, f_result. worst, f_result.f_values, population.get_population());
        // optimizer.gbo();

        VEC_POP x = population.get_population()[0];

        Matrix8x8d modified_dct;
        // Применение преобразования X к DCT-коэффициентам
        population.apply_x_transform(dct_block, x, modified_dct);

        rev_dct_func(image.blocks[i], modified_dct);
        
    }

    int out_width;
    int out_height;

    img_pixels = assemble_image(image, &out_width, &out_height);
    export_image("images/new_lenna.png",img_pixels, 512, 512, 1);


    // Засекаем время окончания выполнения
    auto end = std::chrono::high_resolution_clock::now();

    // Вычисляем разницу между началом и концом в секундах
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start);

    // Выводим время выполнения в секундах
    std::cout << "Time taken by function: " << duration.count() << " seconds" << std::endl;

    return 0; 
}