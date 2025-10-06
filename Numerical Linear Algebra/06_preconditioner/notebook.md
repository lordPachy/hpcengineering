
# Preconditioners with Eigen and LIS

In the following exercise we aim to compare different preconditioning strategies. In Eigen there is a built-in incomplete LU preconditioner that we can provide as an argument for the available linear solver. This ILU preconditioners as two parameters (`fill_in` and `tol`) which can be used to define the accuracy of the approximated LU factorization with respect to the original matrix.\
The general aim is to have
1. $P^{-1}$ or its action should be cheap to compute
2. $\rho(P^{-1}A)$ should be as small as possible (ideally, 1)

Condition 2. arises from the fact that $K(P^{-1}A) \leq \frac{1+\rho(C)}{1-\rho(C)}$\
Some examples of preconditioners are 
 - Jacobi: $P = diag(A)$
 - Sparse Approximate INVerse (SAINV): $P^{-1} = \min_{A\in S\sub \mathbb{R}}(||I - P^{-1}A||^2_F)$, with $||.||_F$ being the Frobenious norm.
 - Incomplete Lower Upper factorization (ILU), where $P = \overline{L} \overline{U}\neq A$  (meaning, approximated); $$\overline{L} = D + L$$ $$\overline{U} = D^{-1}+(D+U)$$
 - General Iterative ILU(P); note that P is the number of iterations allowed
 ```
L0 = U0 = 0
for k = 1...P
  R = A - L0U0
  D = diag(R)
  U0 = uppertri(R)
  L0 = lowertri(R)D^(-1)
end
L = L0 + I
U = U0 + diag(A)
 ```
 Note that if you stop at the first iteration, the sparsity pattern of LU is the same as A, and at each iteration it grows in a squared fashion.
 - ILUT($\epsilon$) (ILU with a threshold): it is performed as actual LU, but entries smaller than $\epsilon$ are discarded
 - Symmetric GS (a sort of ILU(0)) preconditioner is defined, given $A = D + U + L$, as $$P^{-1} = (D + U)^{-1} - (D - U)^{-1} - (D + L)^{-1}=(D+U)^{-1}D(D+L)^{-1}$$
 - Symmetric Sparse OverRelaxation (SSOR), used when A is symmetric $A = D + L + L^T$, is written as above
 - We can parametrize $SSOR(\omega)$, by $$P^{-1} = \frac{\omega}{2-\omega}(\frac{D}{\omega} + L^T)^{-1}D(\frac{D}{\omega} + L)^{-1}$$, for $0 < \omega < 2$
### Exercise 1: Diagonal preconditioner vs. Eigen ILU preconditioner

In the following example, we can observe that the default ILU preconditioner provided by Eigen is very close to a full LU factorization. Indeed, using the `Eigen::IncompleteLUT`, the CG method convergences in just two iterations. 

```
// ... usual include modules
#include "cg.hpp"

int main(int argc, char** argv)
{
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  // Load matrix
  SpMat A;
  Eigen::loadMarket(A, "bcsstm12.mtx");
  // Transform the loaded matrix into a symmetric one
  A = SpMat(A.transpose()) + A;

  double tol = 1.e-14;                 // Convergence tolerance
  int result, maxit = 5000;            // Maximum iterations

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<std::endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<std::endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());

  // Create preconditioners
  Eigen::DiagonalPreconditioner<double> D(A);
  Eigen::IncompleteLUT<double> ILU(A);

  // Hand-made CG with diagonal precond
  x = 0*x;
  result = CG(A, x, b, D, maxit, tol);        // Solve system

  std::cout <<" hand-made CG "<<std::endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: "<<(x-e).norm()<<std::endl;

  // Hand-made CG with ILU precond
  x = 0*x;
  result = CG(A, x, b, ILU, maxit, tol);
  std::cout <<" hand-made CG "<< std::endl;
  cout << "iterations performed: " << maxit << endl;
  cout << "tolerance achieved  : " << tol << endl;
  std::cout << "Error norm: "<<(x-e).norm()<<std::endl;

  return result;
}

```

### Exercise 2: ILU preconditioner in LIS

It is also possible to assess the performances of the ILU preconditioner available in the LIS library and compare it with other preconditioning techniques. Some examples are reported below. 


```
wget https://suitesparse-collection-website.herokuapp.com/MM/HB/bcsstm12.tar.gz
... 

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -maxiter 5000 -tol 1e-12

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-12 -p ilu

mpirun -n 4 ./test1 bcsstm12.mtx 2 sol.txt hist.txt -i bicgstab -tol 1e-12 -p ilu -ilu_fill 2
```
