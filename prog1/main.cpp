#include <iostream>

#include <fstream>

using namespace std;

// Функция для решения системы линейных уравнений методом Гаусса
void solveSystem(double A[4][4], double B[4], double X[4]) {
    // Прямой ход метода Гаусса
    for (int i = 0; i < 3; i++) {
        for (int k = i + 1; k < 4; k++) {
            double factor = A[k][i] / A[i][i];
            B[k] -= factor * B[i];
            for (int j = i; j < 4; j++) {
                A[k][j] -= factor * A[i][j];
            }
        }
    }

    // Обратный ход метода Гаусса
    for (int i = 3; i >= 0; i--) {
        X[i] = B[i] / A[i][i];
        for (int k = i - 1; k >= 0; k--) {
            B[k] -= A[k][i] * X[i];
        }
    }
}

int main() {
    // Задаем координаты точек
    double x_values[] = {M_PI / 4, -M_PI / 4, M_PI, 0};
    double y_values[] = {(M_PI / 4) * (M_PI / 4),
                         (-M_PI / 4) * (-M_PI / 4),
                         M_PI * M_PI,
                        0};

    // Формируем матрицу A и вектор B
    double A[4][4];
    double B[4];
    for (int i = 0; i < 4; ++i) {
        A[i][0] = tan(x_values[i]);
        A[i][1] = (cos(x_values[i]) / 2) + (tan(x_values[i]) / (2 * sqrt(2)));
        A[i][2] = x_values[i] / M_PI;
        A[i][3] = 1 / (2 - sqrt(2));

        B[i] = y_values[i];
    }

    // Решаем систему уравнений
    double X[4];
    solveSystem(A, B, X);

    // Выводим результаты
    cout << "Результаты линейной интерполяции:\n";
    cout << "a: " << X[0] << "\n";
    cout << "b: " << X[1] << "\n";
    cout << "c: " << X[2] << "\n";
    cout << "d: " << X[3] << "\n";

    // Создаем файл данных для построения графика
    ofstream dataFile("data.txt");
    for (double x = -M_PI; x <= M_PI; x += 0.1) {
        double y_interpolation = X[0] * tan(x) + X[1] * ((cos(x) / 2) + (tan(x) / (2 * sqrt(2)))) + X[2] * (x / M_PI) + X[3] * (1 / (2 - sqrt(2)));
        double y_squared = x * x;

        // Проверяем, что y_squared положительно
        if (y_squared > 0) {
            dataFile << x << " " << y_interpolation << " " << y_squared << "\n";
        }
    }
    dataFile.close();

    // Создаем файл данных для точек
    ofstream pointsFile("points.txt");
    for (int i = 0; i < 4; ++i) {
        pointsFile << x_values[i] << " " << y_values[i] << "\n";
    }
    pointsFile.close();

    // Создаем скрипт для gnuplot
    ofstream scriptFile("plot_script.gp");
    scriptFile << "set terminal pngcairo enhanced color size 800,600\n";
    scriptFile << "set output 'interpolation_plot.png'\n";
    scriptFile << "a = " << X[0] << "\n";
    scriptFile << "b = " << X[1] << "\n";
    scriptFile << "c = " << X[2] << "\n";
    scriptFile << "d = " << X[3] << "\n";
    scriptFile << "f(x) = a * tan(x) + b * ((cos(x)/2) + (tan(x)/(2*sqrt(2)))) + c * (x/pi) + d * (1/(2-sqrt(2)))\n";
    scriptFile << "g(x) = x**2\n";
    scriptFile << "set xlabel 'x'\n";
    scriptFile << "set ylabel 'f(x)'\n";
    scriptFile << "set title 'Линейная интерполяция x^2'\n";
    scriptFile << "plot 'data.txt' using 1:2 with lines title 'Интерполяция', \
                             'points.txt' using 1:2 with points title 'Точки' lc rgb 'red', \
                             'data.txt' using 1:3 with lines title 'x^2' lt 2 lc rgb 'blue'\n";
    scriptFile.close();

    // Выполняем gnuplot
    system("gnuplot plot_script.gp");

    return 0;
}
