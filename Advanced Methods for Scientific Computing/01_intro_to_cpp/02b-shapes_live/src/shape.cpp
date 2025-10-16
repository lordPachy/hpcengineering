// It is possible to compile only this file only, then main only, in order to 
// avoid recompiling everything another time
// g++ -C -Wall -I ../include shapes.cpp
// g++ -C -Wall -I ../include main.cpp
// g++ -Wall -I ../include shapes.cc


#include "shape.h"

Circle::Circle(const double r): radius (r) {};
Rectangle::Rectangle(const double a, const double b): height(a), basis(b) {};