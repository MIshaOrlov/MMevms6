//
//  InterpolTriangel.hpp
//  somecode7.5
//
//  Created by Михаил on 31.05.2024.
//

#ifndef InterpolTriangel_hpp
#define InterpolTriangel_hpp

#include "Triangulation.hpp"
#include "Solver.hpp"
#include <stdio.h>
#include <iostream>
#include <cstdlib>
#include <ctime>
#include <cmath>

class Polynomial2D {
public:
    // Конструктор по умолчанию создает нулевой многочлен
    Polynomial2D() {
        for (double& coefficient : coefficients) {
            coefficient = 0;
        }
    }

    // Конструктор, прин&имающий массив коэффициентов
    Polynomial2D(const double coeffs[10]) {
        for (int i = 0; i < 10; ++i) {
            coefficients[i] = coeffs[i];
        }
    }
    
    // Метод для установки коэффициента
    void setCoefficient(int index, double value) {
        if (index >= 0 && index < 10) {
            coefficients[index] = value;
        }
    }
    
    // Метод для получения коэффициента
    double getCoefficient(int index) const {
        if (index >= 0 && index < 10) {
            return coefficients[index];
        }
        return 0; // Возвращаем 0, если индекс некорректный
    }
    
    // Метод для умножения многочлена на число
    Polynomial2D operator*(double multiplier) const {
        Polynomial2D result;
        for (int i = 0; i < 10; ++i) {
            result.setCoefficient(i, coefficients[i] * multiplier);
        }
        return result;
    }
    
    // Метод для деления многочлена на число
    Polynomial2D operator/(double divisor) const {
        Polynomial2D result;
        for (int i = 0; i < 10; ++i) {
            result.setCoefficient(i, coefficients[i] / divisor);
        }
        return result;
    }

    // Оператор присваивания
    Polynomial2D& operator=(const Polynomial2D& other) {
        if (this != &other) {
            for (int i = 0; i < 10; ++i) {
                coefficients[i] = other.coefficients[i];
            }
        }
        return *this;
    }
 
    // Оператор сложения двух многочленов
    Polynomial2D operator+(const Polynomial2D& other) const {
        Polynomial2D result;
        for (int i = 0; i < 10; ++i) {
            result.setCoefficient(i, coefficients[i] + other.coefficients[i]);
        }
        return result;
    }

    // Метод для умножения двух многочленов
    Polynomial2D operator*(const Polynomial2D& other) const {
        double result_coeffs[10] = {0};

        // Перемножаем коэффициенты и складываем результат
        for (int i = 0; i < 10; ++i) {
            for (int j = 0; j < 10; ++j) {
                int deg_x = getDegreeX(i) + getDegreeX(j);
                int deg_y = getDegreeY(i) + getDegreeY(j);
                int index = getIndex(deg_x, deg_y);
                if (index != -1) {
                    result_coeffs[index] += coefficients[i] * other.coefficients[j];
                }
            }
        }

        return Polynomial2D(result_coeffs);
    }

    
    // Метод для вычисления значения многочлена в точке (x, y)
    double evaluate(double x, double y) const {
        double result = 0.0;
        for (int i = 0; i < 10; ++i) {
            result += coefficients[i] * std::pow(x, getDegreeX(i)) * std::pow(y, getDegreeY(i));
        }
        return result;
    }

    
    // Метод для вывода многочлена
    void print() const {
        bool first = true;
        for (int i = 0; i < 10; ++i) {
            if (coefficients[i] != 0) {
                if (!first) {
                    if (coefficients[i] > 0) {
                        std::cout << " + ";
                    } else {
                        std::cout << " - ";
                    }
                } else {
                    if (coefficients[i] < 0) {
                        std::cout << "-";
                    }
                }

                if (std::abs(coefficients[i]) != 1 || (getDegreeX(i) == 0 && getDegreeY(i) == 0)) {
                    std::cout << std::abs(coefficients[i]);
                }

                if (getDegreeX(i) > 0) {
                    std::cout << "x";
                    if (getDegreeX(i) > 1) {
                        std::cout << "^" << getDegreeX(i);
                    }
                }

                if (getDegreeY(i) > 0) {
                    std::cout << "y";
                    if (getDegreeY(i) > 1) {
                        std::cout << "^" << getDegreeY(i);
                    }
                }

                first = false;
            }
        }
        std::cout << std::endl;
    }

private:
    double coefficients[10];

    // Получить степень x для данного индекса коэффициентов
    int getDegreeX(int index) const {
        switch (index) {
            case 0: return 3;
            case 1: return 2;
            case 2: return 1;
            case 3: return 0;
            case 4: return 2;
            case 5: return 1;
            case 6: return 0;
            case 7: return 1;
            case 8: return 0;
            case 9: return 0;
            default: return -1;
        }
    }

    // Получить степень y для данного индекса коэффициентов
    int getDegreeY(int index) const {
        switch (index) {
            case 0: return 0;
            case 1: return 1;
            case 2: return 2;
            case 3: return 3;
            case 4: return 0;
            case 5: return 1;
            case 6: return 2;
            case 7: return 0;
            case 8: return 1;
            case 9: return 0;
            default: return -1;
        }
    }

    // Получить индекс коэффициента по степеням x и y
    int getIndex(int deg_x, int deg_y) const {
        if (deg_x == 3 && deg_y == 0) return 0;
        if (deg_x == 2 && deg_y == 1) return 1;
        if (deg_x == 1 && deg_y == 2) return 2;
        if (deg_x == 0 && deg_y == 3) return 3;
        if (deg_x == 2 && deg_y == 0) return 4;
        if (deg_x == 1 && deg_y == 1) return 5;
        if (deg_x == 0 && deg_y == 2) return 6;
        if (deg_x == 1 && deg_y == 0) return 7;
        if (deg_x == 0 && deg_y == 1) return 8;
        if (deg_x == 0 && deg_y == 0) return 9;
        return -1;
    }
};

double function(double x, double y);

double ValueInterpol(double x, double y,Triangle* triangles, Polynomial2D* phi, Polynomial2D* g,double *alpha, int numTriangles);

void interpolateTriangles(Triangle* triangles, int count,Polynomial2D* g, Polynomial2D* phi, double *alpha);

void triangle_interpol(double x1, double y1, double x2, double y2, double x3, double y3, double* coeff);
    
void psi2(double x1, double y1, double x2, double y2, double x0, double y0, double* coeff);

void approach_triangle(double x1, double y1, double x2, double y2, double x3, double y3, Polynomial2D* g, Polynomial2D* main_phi,int num);

void psi(double u1, double v1, double u2, double v2, Polynomial2D& xyb) ;

bool isInsideTriangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3);

double integrateRectangle(double x0, double x1, double y0, double y1, int numPointsX, int numPointsY, Polynomial2D* main_phi,int j, int k);

double integrateRectangleFunc(double x0, double x1, double y0, double y1, int numPointsX, int numPointsY, Polynomial2D* main_phi,int j) ;

#endif /* InterpolTriangel_hpp */
