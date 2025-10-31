#ifndef DERIVATIVES
#define DERIVATIVES

#include <functional>
#include <type_traits>
namespace DifferenceType
{ 
  struct BACKWARD; 

  struct FORWARD
  {
    using otherType = BACKWARD;
  };

  struct BACKWARD
  {
    using otherType = FORWARD;
  };
}

// RECURSIVE CASE
template <unsigned N, typename F, typename T = double, typename DT = DifferenceType::FORWARD> 
class NthDerivative
{
  // recursion; be careful about the N-1
  public:
    using PreviousDerivative = NthDerivative <N-1,  F, T, typename DT::otherType>;

    // Constructor: we should call the constructors of the previous derivatives
    NthDerivative (const F & f, const T& h_): pDerivative{f, h_}, h(h_) {};

    T operator () (const T& x) const
    { 
      // Doing a type comparison and extracting the value
      if (std::is_same <DifferenceType::FORWARD, DT>::value)
        return (pDerivative(x+h) - pDerivative(x)) / h;
      else
        return (pDerivative(x) - pDerivative(x-h)) / h;
    }

  private:
    PreviousDerivative pDerivative;
    T h;
};

// BASE CASE: SPECIALIZATION
// Note that default template parameters cannot be used in partial specializaition
template <typename F, typename T, typename DT> 
class NthDerivative<1u, F, T, DT>
{
  // recursion; be careful about the N-1
  public:
    // Constructor: we should call the constructors of the previous derivatives
    NthDerivative (const F & f, const T& h_): pDerivative{f}, h(h_) {};

    T operator () (const T& x) const
    { 
      // Doing a type comparison and extracting the value
      if (std::is_same <DifferenceType::FORWARD, DT>::value)
        return (pDerivative(x+h) - pDerivative(x)) / h;
      else
        return (pDerivative(x) - pDerivative(x-h)) / h;
    }

  private:
  // The previous derivative is just the 0-th derivative - the plain function!
    F pDerivative;
    T h;
};

#endif