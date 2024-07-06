//
//  Triangulation.hpp
//  somecode7.5
//
//  Created by Михаил on 31.05.2024.
//

#ifndef Triangulation_hpp
#define Triangulation_hpp


#include <iostream>

struct Point {
    double x, y;
};

struct Triangle {
    Point p1, p2, p3;
};

void triangulateRectangle(Point bottomLeft, Point topRight, int n, Triangle* triangles);

void printTriangles(Triangle* triangles, int count);

int findTriangle(const Point& point, const Triangle* triangles, int numTriangles);

#endif /* Triangulation_hpp */
