//*****************************************************************
// Iterative template routine -- GMRES
//
// GMRES solves the unsymmetric linear system Ax = b using the
// Generalized Minimum Residual method preconditioned with a left preconditioner
//
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
// Input parameters:
//   \tparam Matrix  An Eigen sparse matrix
//   \tparam Vector  An Eigen vector
//   \tparam Preconditioner Something that obeys the preconditioner concept
//   \param A the matrix 
//   \param b the right hand side 
//   \param x the initial guess
//   \param max_iter the maximum number of iterations performed
//   \param m  The restart level
//   \param tol the residual after the final iteration

//*****************************************************************
#include "gmres_util.hpp"
#include <Eigen/Dense>
#include <memory>
namespace LinearAlgebra
{
template <class Matrix, class Vector, class Preconditioner>
int GMRES(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
      int &m, int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Vector::Scalar;
  Real   resid;
  int    i = 0;
  int    j = 1;
  int    k = 0;
  Vector s(m + 1), cs(m + 1), sn(m + 1);
  Vector w(b.size());
  Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic> H =
    Eigen::Matrix<Real, Eigen::Dynamic, Eigen::Dynamic>::Zero(m + 1, m);

  w = b - A * x;
  Vector r = M.solve(w);
  Real   beta = r.norm();
  Real   normb = (beta == 0.0 ? 1.0 : beta); // use preconditioned residual
  resid = r.norm() / normb;
  if(resid <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  std::vector<Vector> v(m + 1, Vector(b.size()));

  while(j <= max_iter)
    {
      v[0] = r * (1.0 / beta); 
      s = Vector::Zero(m + 1);
      s(0) = beta;

      for(i = 0; i < m && j <= max_iter; i++, j++)
        {
          w = M.solve(A * v[i]);
          for(k = 0; k <= i; k++)
            {
              H(k, i) = w.dot(v[k]);
              w -= H(k, i) * v[k];
            }
          H(i + 1, i) = w.norm();
          v[i + 1] = w * (1.0 / H(i + 1, i)); 

          for(k = 0; k < i; k++)
            ApplyPlaneRotation(H(k, i), H(k + 1, i), cs(k), sn(k));

          GeneratePlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
          ApplyPlaneRotation(H(i, i), H(i + 1, i), cs(i), sn(i));
          ApplyPlaneRotation(s(i), s(i + 1), cs(i), sn(i));

          if((resid = std::abs(s(i + 1)) / normb) < tol)
            {
              Update(x, i, H, s, v);
              tol = resid;
              max_iter = j;
              return 0;
            }
        }
      Update(x, m - 1, H, s, v);
      w = b - A * x;
      r = M.solve(w);
      beta = r.norm();
    }

  tol = resid;
  return 1;
}
} // namespace LinearAlgebra
