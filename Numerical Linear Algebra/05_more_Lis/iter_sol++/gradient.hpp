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
int Gradient(const Matrix &A, Vector &x, const Vector &b, const Preconditioner &M,
   int &max_iter, typename Vector::Scalar &tol)
{
  using Real = typename Matrix::Scalar;
  Real   resid;
  Vector r(b.size());
  Vector z(b.size());
  Vector q(b.size());
  Real   alpha;


  Real   normb = b.norm();

  if(normb == 0.0)
    normb = 1;
    

  for(int i = 1; i <= max_iter; i++)
    {   
      r = b - A * x;
      z = M.solve(r);
      //q = A * z;
      alpha = (z.dot(r))/(z.dot(A * z));

      x += alpha * z;
    
      if((resid = r.norm() / normb) <= tol)
        {
          tol = resid;
          std::cout << "Residual: " << resid << std::endl;
          max_iter = i;
          return 0;
        }

      r -= alpha*A * z; 
    }

  tol = resid;
  std::cout << "Residual: " << resid << std::endl;
  return 1;
}
} // namespace LinearAlgebra
