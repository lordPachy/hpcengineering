#include <Derivatives.hpp>
#include <iostream>
#include <cmath>
int main(){

    double h = 0.01;
    std::function<double(const double &) > f = [] (const double & x) {return std::exp(x);};
    std::cout << "The 4th derivative of e^x evaluated in 1 is "
    //<< Nthderivative<...>
    << std::endl;

    return 0;
}