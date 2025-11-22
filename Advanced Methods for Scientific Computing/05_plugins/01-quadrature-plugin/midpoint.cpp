#include "midpoint.hpp"

double integrate(std::function<double (double)> f, double a, double b){
    // Midpoint formula
    return (b-a) * f(0.5*(a+b));
};


#include <cmath>
