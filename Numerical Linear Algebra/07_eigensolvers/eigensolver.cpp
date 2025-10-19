#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    int n = 50;	
    SparseMatrix<double> mat;
    Eigen::loadMarket(mat, "bcsstm12.mtx");

    MatrixXd A;
    A = MatrixXd(mat);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
    if (eigensolver.info() != Eigen::Success) abort();
    std::cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
    // std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
    //           << eigensolver.eigenvectors() << std::endl;
    return 0;    
}

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