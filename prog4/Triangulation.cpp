//
//  Triangulation.cpp
//  somecode7.5
//
//  Created by Михаил on 31.05.2024.
//

#include "Triangulation.hpp"


void triangulateRectangle(Point bottomLeft, Point topRight, int n, Triangle* triangles) {
    double dx = (topRight.x - bottomLeft.x) / n;
    double dy = (topRight.y - bottomLeft.y) / n;
    int count = 0;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            Point p1 = { bottomLeft.x + i * dx, bottomLeft.y + j * dy };
            Point p2 = { bottomLeft.x + (i + 1) * dx, bottomLeft.y + j * dy };
            Point p3 = { bottomLeft.x + i * dx, bottomLeft.y + (j + 1) * dy };
            Point p4 = { bottomLeft.x + (i + 1) * dx, bottomLeft.y + (j + 1) * dy };
            triangles[count].p1 = p1;
            triangles[count].p2 = p2;
            triangles[count].p3 = p4;
            count++;

            triangles[count].p1 = p1;
            triangles[count].p2 = p4;
            triangles[count].p3 = p3;
            count++;
        }
    }
}

void printTriangles(Triangle* triangles, int count) {
    for (int i = 0; i < count; ++i) {
        std::cout << "Triangle " << (i + 1) << ": [("
                  << triangles[i].p1.x << ", " << triangles[i].p1.y << "), ("
                  << triangles[i].p2.x << ", " << triangles[i].p2.y << "), ("
                  << triangles[i].p3.x << ", " << triangles[i].p3.y << ")]\n";
    }
}

int findTriangle(const Point& point, const Triangle* triangles, int numTriangles) {
    for (int i = 0; i < numTriangles; ++i) {
        if((triangles[i].p1.x<= point.x) &&( point.x<=triangles[i].p2.x)){
            if((triangles[i].p1.y<= point.y) &&( point.y<=triangles[i].p3.y)){
                return i;
            }
        }
    }
    return -1; // Точка не находится внутри ни одного треугольника
}
