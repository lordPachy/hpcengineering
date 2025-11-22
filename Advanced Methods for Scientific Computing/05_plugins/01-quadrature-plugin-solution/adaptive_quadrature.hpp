#ifndef HAVE_ADAPTIVE_QUADRATURE_HPP
#define HAVE_ADAPTIVE_QUADRATURE_HPP

#include <functional>

// This is necessary so we can use the function when loaded externally as a single symble.
// It is nice in the case of functions loaded dynamically.
// We must careful with overloading.
extern "C"
{
  double
  integrate(std::function<double(double)>, double, double);
}

#endif
