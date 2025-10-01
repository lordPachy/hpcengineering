//*****************************************************************
// Iterative template routine -- CG
//
// CG solves the symmetric positive definite linear
// system Ax=b using the Conjugate Gradient method.
//
// The return value indicates convergence within max_iter (input)
// iterations (0), or no convergence within max_iter iterations (1).
//
// Upon successful return, output arguments have the following values:
//
//        x  --  approximate solution to Ax = b
// max_iter  --  the number of iterations performed before the
//               tolerance was reached
//      tol  --  the residual after the final iteration
//
//*****************************************************************

namespace LinearAlgebra
{
template <class Matrix, class Vector, class Preconditioner>
int CG(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
   int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Matrix::Scalar;
  Real   resid;
  Vector p(b.size());
  Vector z(b.size());
  Vector q(b.size());
  Real   alpha, beta, rho;
  Real   rho_1(0.0);

  Real   normb = b.norm();
  Vector r = b - A * x;

  if(normb == 0.0)
    normb = 1;

  if((resid = r.norm() / normb) <= tol)
    {
      tol = resid;
      max_iter = 0;
      return 0;
    }

  for(int i = 1; i <= max_iter; i++)
    {
      z = M.solve(r);
      rho = r.dot(z);

      if(i == 1)
        p = z;
      else
        {
          beta = rho / rho_1;
          p = z + beta * p;
        }

      q = A * p;
      alpha = rho / p.dot(q);

      x += alpha * p;
      r -= alpha * q;

      if((resid = r.norm() / normb) <= tol)
        {
          tol = resid;
          max_iter = i;
          return 0;
        }

      rho_1 = rho;
    }

  tol = resid;
  return 1;
}
} // namespace LinearAlgebra
