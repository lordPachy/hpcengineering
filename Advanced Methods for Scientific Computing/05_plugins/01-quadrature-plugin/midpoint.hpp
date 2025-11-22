#ifndef HAVE_MIDPOINT_HPP
#define HAVE_MIDPOINT_HPP

#include <functional>

double integrate(std::function<double (double)> f, double a, double b);

#endif
