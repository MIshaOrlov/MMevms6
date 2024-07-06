//
//  InterpolTriangel.cpp
//  somecode7.5
//
//  Created by Михаил on 31.05.2024.
//

#include "InterpolTriangel.hpp"


double function(double x, double y){
    //(x*x*y)*sin(x*y);
    return std::sin(M_PI * x) * std::cos(M_PI * y);
}


double ValueInterpol(double x, double y,Triangle* triangles, Polynomial2D* phi, Polynomial2D* g,double *alpha, int numTriangles){
    Point p = {x,y};
    double coeffs[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int temp = findTriangle(p, triangles,  numTriangles);
    
    if(temp < 0){
        std::cout<<x<<" "<<y<<std::endl;
        return -1;
    } else {
        return g[temp].evaluate(x, y);
    }
}


void interpolateTriangles(Triangle* triangles, int count,Polynomial2D* g,Polynomial2D* phi, double *alpha) {
    double A[100] = {0.0};
    //double * coeff = new double[10];
    double B[10]= {0.0};
    double x[10];
    for (int i = 0; i < count; i++) {
        // g многочлен на i треугольнике вычисляется получается
        if(i % 2 == 0){
            approach_triangle(triangles[i].p1.x, triangles[i].p1.y, triangles[i].p2.x, triangles[i].p2.y, triangles[i].p3.x, triangles[i].p3.y,g,phi,i);
            //triangle_interpol(triangles[i].p1.x, triangles[i].p1.y, triangles[i].p2.x, triangles[i].p2.y, triangles[i].p3.x, triangles[i].p3.y, coeff);
        } else {
            approach_triangle(triangles[i].p1.x, triangles[i].p1.y, triangles[i].p3.x, triangles[i].p3.y, triangles[i].p2.x, triangles[i].p2.y,g,phi,i);
            //triangle_interpol(triangles[i].p1.x, triangles[i].p1.y, triangles[i].p3.x, triangles[i].p3.y, triangles[i].p2.x, triangles[i].p2.y, coeff);
        }
        
        for(int j = 0;j < 10;j++){
            for(int k = j; k < 10;k++){
                
                A[j*10 + k] = integrateRectangle(triangles[i].p1.x, triangles[i].p2.x,  triangles[i].p1.y, triangles[i].p3.y, 100, 100, phi,j,k );
                A[j + k*10] = A[j*10 + k];
            }
        }
        for(int j = 0; j < 10; j++){
            B[j] =  integrateRectangleFunc(triangles[i].p1.x, triangles[i].p2.x,  triangles[i].p1.y, triangles[i].p3.y, 100, 100, phi,j );
        }
        
       solver(10, A, B, x);
        for(int j = 0;j < 10;j++){
            alpha[i+j] = x[j];
        }
    }
    
}

// Функция для вычисления значения в точке x4,y4
void approach_triangle(double x1, double y1, double x2, double y2, double x3, double y3, Polynomial2D* g, Polynomial2D* main_phi,int num){
    // Вычисляем координаты точек, делящих стороны треугольника на три равные части
    double x[10], y[10],temp;
    x[0] = x1; y[0] = y1;
    x[1] = x2; y[1] = y2;
    x[2] = x3; y[2] = y3;
    x[3] = (2 * x1 + x2) / 3.0; y[3] = (2 * y1 + y2) / 3.0;
    x[4] = (x1 + 2 * x2) / 3.0; y[4] = (y1 + 2 * y2) / 3.0;
    x[5] = (2 * x2 + x3) / 3.0; y[5] = (2 * y2 + y3) / 3.0;
    x[6] = (x2 + 2 * x3) / 3.0; y[6] = (y2 + 2 * y3) / 3.0;
    x[7] = (2 * x1 + x3) / 3.0; y[7] = (2 * y1 + y3) / 3.0;
    x[8] = (x1 + 2 * x3) / 3.0; y[8] = (y1 + 2 * y3) / 3.0;
    x[9] = (x1 + x2 + x3) / 3.0; y[9] = (y1 + y2 + y3) / 3.0;
    double coeff[3];
    // Вычисляем пси
    double coeffs[10] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    Polynomial2D zero(coeffs) ;
    Polynomial2D psi_[10];
    for(int i = 0; i< 10; i++){
        psi_[i] = zero;
    }
    
     psi(x[1],y[1],x[2],y[2],psi_[0]);
     psi(x[0],y[0],x[2],y[2],psi_[1]);
     psi(x[0],y[0],x[1],y[1],psi_[2]);
     psi(x[3],y[3],x[7],y[7],psi_[3]);
     psi(x[4],y[4],x[8],y[8],psi_[4]);
     psi(x[4],y[4],x[5],y[5],psi_[5]);
     psi(x[3],y[3],x[6],y[6],psi_[6]);
     psi(x[6],y[6],x[8],y[8],psi_[7]);
     psi(x[5],y[5],x[7],y[7],psi_[8]);
    
    Polynomial2D phi[10];
    for(int i = 0; i< 10; i++){
        phi[i] = zero;
    }
    phi[0] = psi_[0]*psi_[3]*psi_[4];
    phi[1] = psi_[1]*psi_[5]*psi_[6];
    phi[2] = psi_[2]*psi_[7]*psi_[8];
    phi[3] = psi_[0]*psi_[1]*psi_[4];
    //std::cout<<"Checkout"<<std::endl;
    //psi_[0].print();
    //std::cout<<"Check " << psi_[0].evaluate(x[2], y[2]) <<std::endl;
   // psi_[1].print();
   // psi_[4].print();
    //phi[3].print();
    //std::cout<<"Checkout " <<phi[3].evaluate(x[3], y[3]) <<std::endl;
    phi[4] = psi_[0]*psi_[1]*psi_[6]; phi[5] = psi_[1]*psi_[2]*psi_[6];
    phi[6] = psi_[1]*psi_[2]*psi_[8]; phi[7] = psi_[0]*psi_[2]*psi_[4];
    phi[8] = psi_[0]*psi_[2]*psi_[8]; phi[9] = psi_[0]*psi_[1]*psi_[2];
    for(int i = 0; i < 10; i++){
        double temp = (phi[i].evaluate(x[i], y[i]));
        phi[i] = phi[i]/temp;
        //phi[i].print();
        //std::cout<<"Check " <<phi[i].evaluate(x[3], y[3]) <<std::endl;
        main_phi[num + i] = phi[i];
        g[num] = g[num] + phi[i]*function(x[i], y[i]);
    }
    //std::cout<<"Checkout end"<<std::endl;
}


/*L(u1,v1), (u2,v2) (x, y)*/
void psi(double u1, double v1, double u2, double v2,Polynomial2D& xyb){
    
    xyb.setCoefficient(7,v2 - v1);
    xyb.setCoefficient(8,-(u2 - u1));
    xyb.setCoefficient(9,v1*(u2 - u1) - u1*(v2 - v1));
}


// Определяет, находится ли точка (x, y) внутри треугольника
bool isInsideTriangle(double x, double y, double x1, double y1, double x2, double y2, double x3, double y3) {
    double denominator = ((y2 - y3) * (x1 - x3) + (x3 - x2) * (y1 - y3));
    double a = ((y2 - y3) * (x - x3) + (x3 - x2) * (y - y3)) / denominator;
    double b = ((y3 - y1) * (x - x3) + (x1 - x3) * (y - y3)) / denominator;
    double c = 1 - a - b;
    return 0 <= a && a <= 1 && 0 <= b && b <= 1 && 0 <= c && c <= 1;
}

// Функция для вычисления двумерного интеграла методом прямоугольников
double integrateRectangle(double x0, double x1, double y0, double y1, int numPointsX, int numPointsY, Polynomial2D* main_phi,int j, int k) {
    double deltaX = (x1 - x0) / numPointsX;
    double deltaY = (y1 - y0) / numPointsY;
    double integral = 0.0;

    for (int i = 0; i < numPointsX; ++i) {
        for (int j = 0; j < numPointsY; ++j) {
            double xLeft = x0 + i * deltaX;
            double xRight = x0 + (i + 1) * deltaX;
            double yBottom = y0 + j * deltaY;
            double yTop = y0 + (j + 1) * deltaY;
            double fValue = main_phi[i].evaluate((xLeft + xRight) / 2.0, (yBottom + yTop) / 2.0) * main_phi[k].evaluate((xLeft + xRight) / 2.0, (yBottom + yTop) / 2.0);
            double area = (xRight - xLeft) * (yTop - yBottom);

            integral += fValue * area;
        }
    }

    return integral;
}

// Функция для вычисления двумерного интеграла методом прямоугольников
double integrateRectangleFunc(double x0, double x1, double y0, double y1, int numPointsX, int numPointsY, Polynomial2D* main_phi,int j) {
    double deltaX = (x1 - x0) / numPointsX;
    double deltaY = (y1 - y0) / numPointsY;
    double integral = 0.0;

    for (int i = 0; i < numPointsX; ++i) {
        for (int j = 0; j < numPointsY; ++j) {
            double xLeft = x0 + i * deltaX;
            double xRight = x0 + (i + 1) * deltaX;
            double yBottom = y0 + j * deltaY;
            double yTop = y0 + (j + 1) * deltaY;
            double fValue = main_phi[i].evaluate((xLeft + xRight) / 2.0, (yBottom + yTop) / 2.0) * function((xLeft + xRight) / 2.0, (yBottom + yTop) / 2.0);
            double area = (xRight - xLeft) * (yTop - yBottom);

            integral += fValue * area;
        }
    }

    return integral;
}
