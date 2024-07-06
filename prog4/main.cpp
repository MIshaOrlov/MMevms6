/*
 Приближение функции, заданной в прямоугольнике, конечными элементами степени 3 методом наи- меньших квадратов.
 a1*x^3 + a2*x^2*y + a3*x*y^2 + a4*y^3 + a5*x^2 + a6*x*y + a7*y^2 + a8*x + a9*y + a10

 */


#include <iostream>
#include <cmath>
#include "Triangulation.hpp"
#include "InterpolTriangel.hpp"
#include <cstdlib>
#include <fstream>
#include <ctime>

using namespace std;

double f(double x, double y){
    return std::sin(M_PI * x) * std::cos(M_PI * y);
}



int main(int argc, char* argv[]) {
    if (argc != 6) { // Нужно 5 аргументов (n и 4 координаты двух точек)
            std::cout << "Usage: " << argv[0] << " x1 y1 x2 y2 n" << std::endl;
            return 1;
        }
    int n;
    double x1,y1,x2,y2;
    try {
        x1 = std::stod(argv[1]);
        y1 = std::stod(argv[2]);
        x2 = std::stod(argv[3]);
        y2 = std::stod(argv[4]);
        n = std::stoi(argv[5]);
    } catch (const std::invalid_argument& e) {
        std::cout << "Invalid argument: please provide numeric values for coordinates." << std::endl;
        return 1;
    } catch (const std::out_of_range& e) {
        std::cout << "Argument out of range: please provide realistic numeric values for coordinates." << std::endl;
        return 1;
    }
    std::cout << "The value of n is: " << n << std::endl;
    std::cout << "The coordinates of the rectangle are: (" << x1 << ", " << y1 << ") and (" << x2 << ", " << y2 << ")" << std::endl;

    Point bottomLeft = {x1, y1};
    Point topRight = {x2, y2};
    
    
    
    std::clock_t start = std::clock();
    Triangle* triangles = new Triangle[2 * n * n];
    triangulateRectangle( bottomLeft, topRight, n, triangles);
    printTriangles(triangles, 2*n*n);
    std::clock_t end = std::clock();
    double duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Время выполнения триангуляции: " << duration << " секунд" << std::endl;

    //printTriangles(triangles, 2 * n * n);
    
    start = std::clock();
    Polynomial2D* g  = new Polynomial2D[2* n * n];
    Polynomial2D* phi  = new Polynomial2D[2* n * n * 10];
    double* apha = new double[2* n * n * 10];
    interpolateTriangles(triangles , 2 * n * n , g, phi,apha );
    end = std::clock();
    duration = double(end - start) / CLOCKS_PER_SEC;
    std::cout << "Время выполнения поиска коэффициентов: " << duration << " секунд" << std::endl;
    
    /*for(int i = 0; i < 2* n * n * 10; i++ ){
        for(int i = j; j < 2* n * n * 10; j++ ){
            phi[i].evaluate(x[j], y[j]);
        }
    }*/
    
    // Генерация данных для функции f(x, y)
    std::ofstream function_file("function.dat");
    double x_min = 0.0, x_max = 3.0;
    double y_min = 0.0, y_max = 3.0;
    int grid_size = 100; // Количество точек в сетке по каждой оси
 // Количество точек в сетке по каждой оси
    for (int i = 0; i <= grid_size; ++i) {
        for (int j = 0; j <= grid_size; ++j) {
            double x = x_min + i * (x_max - x_min) / grid_size;
            double y = y_min + j * (y_max - y_min) / grid_size;
            double z = f(x, y);
            function_file << x << " " << y << " " << z << "\n";
        }
        function_file << "\n"; // разделитель строк для gnuplot
    }
    function_file.close();
    
    // Генерация данных для функции f(x, y)
    std::ofstream interpol_file("interpol.dat");
    for (int i = 0; i <= grid_size; ++i) {
        for (int j = 0; j <= grid_size; ++j) {
            double x = x_min + i * (x_max - x_min) / grid_size;
            double y = y_min + j * (y_max - y_min) / grid_size;
            double z = ValueInterpol(x, y, triangles, phi,g, apha, 2*n*n);
            interpol_file << x << " " << y << " " << z << "\n";
        }
        interpol_file << "\n"; // разделитель строк для gnuplot
    }
    interpol_file.close();

    // Создание скрипта для Gnuplot
    std::ofstream gnuplot_script("triangulation.gp");
    gnuplot_script << "set title \"Triangulation of a Rectangle and Function f(x, y)\"\n";
    gnuplot_script << "set xrange [0:3]\n";
    gnuplot_script << "set yrange [0:3]\n";
    gnuplot_script << "set zrange [-2:2]\n";
    gnuplot_script << "set size ratio -1\n";
    gnuplot_script << "unset key\n";
    gnuplot_script << "set dgrid3d 50,50\n";
    gnuplot_script << "set hidden3d\n";
    //gnuplot_script << "splot 'points.dat' using 1:2:(0) with points pt 7 lc rgb \"red\" title 'Points', \\\n";
    //gnuplot_script << "      'triangles.dat' using 1:2:(0) with lines lc rgb \"blue\" title 'Triangles', \\\n";
    gnuplot_script << "splot   'function.dat' using 1:2:3 with lines lc rgb \"green\" title 'Function f(x, y)', \\\n";
    gnuplot_script << "        'interpol.dat' using 1:2:3 with lines lc rgb \"blue\" title 'Interpol Function f(x, y)'\n";
    gnuplot_script << "pause -1 \n";
    gnuplot_script.close();
    
    std::cout << "Script file 'triangulation.gp' has been created.\n";
    std::cout << "You can now run 'gnuplot triangulation.gp' to generate the plot.\n";
    
    
    delete[] triangles;
    delete [] apha;
    delete [] phi;
    return 0;
}


