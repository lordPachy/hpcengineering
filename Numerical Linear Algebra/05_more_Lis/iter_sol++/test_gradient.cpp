#include <cstdlib>                      // System includes
#include <iostream>
#include <Eigen/SparseCore>
#include <Eigen/IterativeLinearSolvers>
#include <chrono>

using std::endl;
using std::cout;
#include "gradient.hpp"                       // Include the template c
int main(int argc, char** argv){
  using namespace LinearAlgebra;
  // Some useful alias
  using SpMat=Eigen::SparseMatrix<double>;
  using SpVec=Eigen::VectorXd;

  int n = 1000;
  SpMat A(n,n);                       // define matrix
  A.reserve(2998);
  for (int i=0; i<n; i++) {
      A.coeffRef(i, i) = 2.0*(i+1);
      if(i>0) A.coeffRef(i, i-1) = -i;
      if(i<n-1) A.coeffRef(i, i+1) = -(i+1);
  }
  double tol = 1.e-10;                // Convergence tolerance
  int result, maxit = 50000;           // Maximum iterations

  std::cout<<"Matrix size:"<<A.rows()<<"X"<<A.cols()<<endl;
  std::cout<<"Non zero entries:"<<A.nonZeros()<<endl;

  SpMat B = SpMat(A.transpose()) - A;  // Check symmetry
  std::cout<<"Norm of A-A.t: "<<B.norm()<<endl;

  // Create Rhs b
  SpVec e = SpVec::Ones(A.rows());
  SpVec b = A*e;
  SpVec x(A.rows());
  Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner
  // It is important to use this type as it allows the method .solve()
  // First with eigen CG
  Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
  cg.setMaxIterations(maxit);
  cg.setTolerance(tol);
  cg.compute(A);
  x = cg.solve(b);
  std::cout << "Eigen standard CG method"<< endl;
  std::cout << "#iterations:     " << cg.iterations() << endl;
  std::cout << "estimated error: " << cg.error()      << endl;
  std::cout << "effective error: "<<(x-e).norm()<< endl;

  // Now with hand-made gradient
  x=0*x;
  std::cout <<" Hand-made gradient method "<< endl;

  auto start = std::chrono::high_resolution_clock::now();
  result = Gradient(A, x, b, D, maxit, tol);        // Solve system
  auto stop = std::chrono::high_resolution_clock::now();

  auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(stop - start);
  cout << "Elapsed time: " << duration.count() << endl;
  cout << "Gradient flag = " << result << endl;
  cout << "iterations performed: " << maxit << endl;

  return result;
}