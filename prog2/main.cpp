#include <iostream>
#include <fstream>
#include <vector>
#include <cmath>

using namespace std;


const double goldenRatio = (1 + sqrt(5)) / 2;

double funct(double x){
    return x*x*x*x*x*x + x*x*sin(x);
}

double derevative_function(double x){
    return 6*x*x*x*x*x + 2*x*sin(x) + x*x*cos(x);
}

double secondDerivative(double x) {
    return 30 * x*x*x*x + 2 * sin(x) + 4 * x * cos(x) + 2  *x*x * cos(x) - x*x * sin(x);
}

double fourth_derivative(double x){
    return  360*x*x - 12*sin(x) - 8*x*cos(x)+sin(x)*x*x;
}

// Функция для кубической интерполяции Эрмита
double cubicHermiteInterpolation(double x, const double* x_values, const double* y_values, const double* dy_values, int size) {
    // Находим интервал, в который попадает x
    int interval = 0;
    while (interval < size - 1 && x > x_values[interval + 1])
        interval++;

    // Вычисляем локальный параметр t в интервале
    double t = (x - x_values[interval]) / (x_values[interval + 1] - x_values[interval]);
    double t2 = t * t;
    double t3 = t2 * t;

    // Вычисляем значения полинома Эрмита
    double h00 = 2 * t3 - 3 * t2 + 1;
    double h10 = t3 - 2 * t2 + t;
    double h01 = -2 * t3 + 3 * t2;
    double h11 = t3 - t2;

    // Вычисляем значение функции в точке x
    double interpolated_value = h00 * y_values[interval] +
                                     h10 * dy_values[interval] * (x_values[interval + 1] - x_values[interval]) +
                                     h01 * y_values[interval + 1] +
                                     h11 * dy_values[interval + 1] * (x_values[interval + 1] - x_values[interval]);

    return interpolated_value;
}

// Функция для нахождения максимума на отрезке [a, b]
double findMax(double a, double b, double epsilon,  const double* x_values, const double* y_values, const double* dy_values, int size) {
    double x1, x2;
    
    while (fabs(b - a) > epsilon) {
           x1 = b - (b - a) / goldenRatio;
           x2 = a + (b - a) / goldenRatio;
           if (fabs(funct(x1) - cubicHermiteInterpolation(x1, x_values, y_values, dy_values,size)) <= fabs(funct(x2) - cubicHermiteInterpolation(x2, x_values, y_values, dy_values,size))) // Условие для поиска максимума
               a = x1;
           else
             b = x2;
    } // Выполняем, пока не достигнем заданной точности
    return (a + b)/2;
}

double  findMaxDerevative(double a, double b, double epsilon) {
    double x1, x2;
    
    while (fabs(b - a) > epsilon) {
           x1 = b - (b - a) / goldenRatio;
           x2 = a + (b - a) / goldenRatio;
           if (fourth_derivative(x1) <= fourth_derivative(x2)) // Условие для поиска максимума
               a = x1;
           else
             b = x2;
    } // Выполняем, пока не достигнем заданной точности
    return (a + b)/2;
}




int main() {
    // Входные данные
    const int size = 5;
    double* x_values = new double[size];
    double* y_values = new double[size];
    double* dy_values = new double[size];
    double max_diff = 0.0;
    double temp = 0.0;
    double x_cr = 0.0;
    double y_ct = 0.0,y_cr = 0.0;
    double fourth_derivative_max = 0.0;

       // Заполняем массивы значениями
    for (int i = 0; i < size; ++i) {
        x_values[i] = i+1;
        y_values[i] = funct(i+1);
        dy_values[i] = derevative_function(i+1);
    }
    //Вычисляем верхнюю оценку для ошибки интерполяции
    double h_max = 0.0;
    for (size_t i = 0; i < size - 1; ++i) {
        double h = x_values[i + 1] - x_values[i];
        if (h > h_max) {
            h_max = h;
        }
    }
    

    // Создаем файл данных для построения графика
    std::ofstream dataFile("data.txt");
    for (double x = 1; x <= 5; x += 0.01) {
        double y_interpolation = cubicHermiteInterpolation(x, x_values, y_values, dy_values,size);
        double y_squared = funct(x);
        dataFile << x << " " << y_interpolation << " " << y_squared << "\n";
    }
    dataFile.close();
    
    
    
    //найдем максимум разности
    x_cr = findMax(4,5,0.0001, x_values, y_values, dy_values,size);
    y_ct = funct(x_cr);
    y_cr = cubicHermiteInterpolation(x_cr, x_values, y_values, dy_values,size);
    max_diff = fabs(y_ct - y_cr);
    
    fourth_derivative_max = fourth_derivative(findMaxDerevative(1,5,0.0001));
    //cout<< fourth_derivative_max<< endl;
    
    double upper_bound = (1.0 / (16*24)) * h_max*h_max*h_max*h_max * fourth_derivative_max;
    
    
    // Создаем файл данных для точек
    std::ofstream pointsFile("points.txt");
    for (size_t i = 0; i < size; ++i) {
        pointsFile << x_values[i] << " " << y_values[i] << "\n";
    }
    pointsFile << x_cr << " " << y_ct << "\n";
    pointsFile << x_cr << " " << y_cr << "\n";
    pointsFile.close();
    
   
    // Создаем скрипт для gnuplot
    std::ofstream scriptFile("plot_script.gp");
    scriptFile << "set terminal pngcairo enhanced color size 800,600\n";
    scriptFile << "set output 'interpolation_plot.png'\n";
    scriptFile << "set xlabel 'x'\n";
    scriptFile << "set ylabel 'f(x)'\n";
    scriptFile << "set title 'Интерполяция кубическими сплайнами Эрмита'\n";
    scriptFile << "set xrange [0:6]\n";
    scriptFile << "plot 'data.txt' using 1:2 with lines lw 2 lc rgb 'blue' title 'Интерполяция', "
                      "'points.txt' using 1:2 with points title 'Точки', "
                      "'data.txt' using 1:3 with lines lw 1 lc rgb 'red' title 'x^6 + x^2 * sin(x)'\n";
    scriptFile.close();

    // Выполняем gnuplot
    system("gnuplot plot_script.gp");
    
    cout<<"Максимальная разница на графике : "<< max_diff<<endl;
    cout<<"Точка на графике которая проблема : "<< x_cr <<" ориг : "<< y_ct << " интерпол : " << y_cr << endl ;
    cout<<"Значение верхней оценки для соответствующего метода интерполяции : "<<upper_bound<<endl;

    // Создаем скрипт для gnuplot
    std::ofstream diffScriptFile("diff_plot_script.gp");
    diffScriptFile << "set terminal pngcairo enhanced color size 800,600\n";
    diffScriptFile << "set output 'difference_plot.png'\n";
    diffScriptFile << "set xlabel 'x'\n";
    diffScriptFile << "set ylabel 'f(x) - Interpolation'\n";
    diffScriptFile << "set title 'Разница между f(x) и интерполяцией'\n";
    diffScriptFile << "set xrange [0:6]\n";
    diffScriptFile << "plot 'difference_data.txt' using 1:2 with lines lw 2 lc rgb 'green' title 'Разница'\n";
    diffScriptFile.close();

    // Создаем файл данных для разницы
    std::ofstream diffDataFile("difference_data.txt");
    for (double x = 1; x <= 5; x += 0.01) {
        double y_difference = funct(x) - cubicHermiteInterpolation(x, x_values, y_values, dy_values,size);
        diffDataFile << x << " " << y_difference << "\n";
    }
    diffDataFile.close();

    // Выполняем gnuplot для построения графика разницы
    system("gnuplot diff_plot_script.gp");
    
    return 0;
}
