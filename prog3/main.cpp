//
//  main.cpp
//  somecode7.6
//
//  Created by Михаил on 29.05.2024.
//

/*
 Интерполяция параболическими сплайнами с определением недостающих граничных условий при по- мощи дополнительного узла в приграничных узлах.
 */

#include <iostream>
#include <fstream>


void printArray(const double* arr, int size) {
    std::cout << "Массив: ";
    for (int i = 0; i < size; ++i) {
        std::cout << arr[i] << " ";
    }
    std::cout << std::endl;
}

void generate_unique_points(double a, double b, int n, double* points) {

    // Генерация точек
    int count = 0;
    while (count < n) {
        double point = a + count*(b-a)/n;
        points[count++] = point;
    }
}

double function(double x){
    return x*x;
}

/*
 
 c[0] /= a1[0];
 for (int i = 1; i < n - 1; ++i) {
     a1[i] -= a3[i - 1] * c[i - 1];
     c[i] /= a1[i];
 }
 a1[n - 1] -= a3[n - 2] * c[n - 2];

 x[0] = b[0] / a1[0];
 for (int i = 1; i < n; ++i)
     x[i] = (b[i] - a3[i - 1] * x[i - 1]) / a1[i];

 for (int i = n - 2; i >= 0; --i)
     x[i] -= c[i] * x[i + 1];

 delete[] c;
 delete[] d;

 */
void solveTridiagonal(int n, double *a, double *c, double *d, double *b, double *x)
{
    int i;

    c[0] /= a[0];
    for (i = 1; i < n - 1; i++)
    {
        a[i] -= d[i - 1] * c[i - 1];
        c[i] /= a[i];
    }
    a[n - 1] -= d[n - 2] * c[n - 2];

    x[0] = b[0] / a[0];
    for (i = 1; i < n; i++)
        x[i] = (b[i] - d[i - 1] * x[i - 1]) / a[i];

    for (i = n - 2; i >= 0; i--)
        x[i] -= c[i] * x[i + 1];
}

int  interpol(int n, double * ksi, double*x, double* y, double* c1,double* c2, double* c3);

double Pfunction(double* c1,double* c2, double* c3, double* ksi, double x,int n){
    for(int i = 0; i < n; i++){
        if(x < ksi[i+1]){
            return c1[i]+c2[i]*(x-ksi[i]) + c3[i]*(x-ksi[i])*(x-ksi[i]);
        }
    }
    return -1;
}

int main() {
    int n,a,b;
    std::cin>>n>>a>>b;
    double* x = new double[n];
    double* y = new double[n];
    double* ksi = new double[n+1];
    double* c1 = new double[n];
    double* c2 = new double[n];
    double* c3 = new double[n];
    
    generate_unique_points(a,b,n,x);
    
    interpol(n, ksi, x, y, c1, c2, c3);
    
    std::cout<<1<<std::endl;
    // Создаем файл данных для построения графика
    std::ofstream dataFile("data.txt");
    for (double t = a  ; t <= b; t += 0.1) {

        dataFile << t << " " << Pfunction(c1,c2,c3,ksi,t,n)<< " " << function(t)<< "\n";
    }
    dataFile.close();

    
    // Создаем скрипт для gnuplot
    std::ofstream scriptFile("plot_script.gp");
    scriptFile << "set terminal pngcairo enhanced color size 800,600\n";
    scriptFile << "set output 'interpolation_plot.png'\n";
    scriptFile << "set xlabel 'x'\n";
    scriptFile << "set ylabel 'f(x)'\n";
    scriptFile << "set xrange [-10:10]\n";
    scriptFile << "set yrange [-10:10]\n";
    scriptFile << "set title 'Линейная интерполяция x^2'\n";
    scriptFile << "plot 'data.txt' using 1:2 with lines title 'Интерполяция' lw 5, \
                             'data.txt' using 1:3 with lines title 'function' lt 2 lc rgb 'blue'\n";
    scriptFile.close();

    // Выполняем gnuplot
    system("gnuplot plot_script.gp");

    // Уведомление об успешном завершении
    std::cout << "Графики функций сохранены в файле functions_plot.png" << std::endl;

    
    
    delete[] x;
    delete[] y;
    delete[] ksi;
    delete[] c1;
    delete[] c2;
    delete[] c3;
    return 0;
}

int  interpol(int n, double * ksi, double*x, double* y, double* c1,double* c2, double* c3){
    double* koef_1 = new double[n];
    double* koef_2 = new double[n+1];
    double* koef_3= new double[n];
    double* koef_4= new double[n];
    double* v = new double[n+1];
    
    for (int i = 1; i < n; i++)
        ksi[i] = (x[i - 1] + x[i])/2;
    
    ksi[0] = x[0] - (x[1])/2;
    ksi[n] = x[n-1] + (x[n-2])/2;
    
    for (int i = 1; i < n; i++)
    {
        koef_1[i-1] = (1.0/(x[i-1]-ksi[i-1]) - 1.0/(ksi[i]-ksi[i-1]));
        koef_2[i] = (1.0/(ksi[i] - x[i - 1]) + 1.0/(ksi[i] - ksi[i - 1]) + 1.0/(x[i] - ksi[i]) + 1.0/(ksi[i + 1] - ksi[i]));
        koef_3[i] =  (1.0/(ksi[i + 1] - x[i]) - 1.0/(ksi[i + 1] - ksi[i]));
        koef_4[i] = koef_1[i-1]*function(x[i-1]) +(1.0/(x[i] - ksi[i]) + 1.0/(ksi[i + 1] - x[i]))*function(x[i]);
    }
    
    // Заполняем граничные условия
    
    koef_2[0] = 1.0;
    koef_3[0] = 0.0;
    koef_1[n-1] = 0.0;
    koef_2[n] = 1.0;
    koef_4[0] = function(ksi[0]);
    koef_4[n] = function(ksi[n]);
    
    //printArray(koef_1,n);
    //printArray(koef_2,n+1);
    //printArray(koef_3,n);
    //printArray(koef_4,n);
    solveTridiagonal(n + 1, koef_2, koef_3, koef_1, koef_4, v);
    for(int i = 0;i < n; i++){
        /*c1[i] = v[i];
        c2[i] = ((function(x[i]) - v[i])/(x[i] - ksi[i])) - ((x[i] - ksi[i])/(ksi[i+1] - ksi[i]))*((v[i+1]-function(x[i]))/(ksi[i+1] - x[i]) - (function(x[i]) - v[i])/(x[i] - ksi[i]));
        c3[i] =(1/(ksi[i+1] - ksi[i]))*((v[i+1]-function(x[i]))/(ksi[i+1] - x[i]) - (function(x[i]) - v[i])/(x[i] - ksi[i]));*/
        c1[i] = v[i];
        double tmp1 = ((v[i + 1] - function(x[i])) / (ksi[i + 1] - x[i]) - (function(x[i]) - v[i]) / (x[i] - ksi[i])) / (ksi[i + 1] - ksi[i]);
        c2[i] = (function(x[i]) - v[i]) / (x[i] - ksi[i]) - (x[i] - ksi[i]) * tmp1;
        c3[i] = tmp1;
        
    }
    //printArray(c1,n);
    //printArray(c2,n);
    //printArray(c3,n);
    
    delete [] koef_1;
    delete [] koef_2;
    delete [] koef_3;
    delete [] koef_4;
    return 1;
}

