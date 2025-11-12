#include <iostream>
#include <unsupported/Eigen/SparseExtra>
#include <Eigen/Sparse>
#include "grad.hpp"
#include "cg.hpp"

using namespace Eigen;
using namespace std;

using SpMat = SparseMatrix<double>;
using SpVec = VectorXd;

int main(){
    SpMat A;
    loadMarket(A, "Aex1.mtx");

    // Checking symmetry and computing norm of A
    SpMat A_t = A.transpose();
    cout << "Norm of A - A_T: " << (A - A_t).norm() << endl;
    cout << "Norm of A: " << A.norm() << endl;

    // Computing the product
    SpVec x_e(A.cols());
    for (int i = 0; i < A.cols(); i++){
        if (i % 2 == 0){
            x_e.coeffRef(i) = 1;
        } else {
            x_e.coeffRef(i) = -1;
        }
    }

    SpVec b = A*x_e;
    cout << "Norm of b: " << b.norm() << endl;

    // Solving the system
    double tol = 1.e-11;                  // Convergence tolerance
    int result, maxit = 20000;           // Maximum iterations
    SpVec x(A.rows());
    Eigen::DiagonalPreconditioner<double> D(A);

    // Gradient method solution
    x=0*x;
    result = LinearAlgebra::GRAD(A, x, b, D, maxit, tol);        // Solve system

    std::cout <<"GRADIENT METHOD"<< endl;
    cout << "Gradient flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    std::cout << "Relative error norm: "<<((x-x_e).norm())/(x_e.norm())<< endl;

    // Conjugate Gradient Solution
    x=0*x;
    result = LinearAlgebra::CG(A, x, b, D, maxit, tol);        // Solve system

    std::cout <<"CONJUGATE GRADIENT"<< endl;
    cout << "CG flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    std::cout << "Relative error norm: "<<(x-x_e).norm()/(x_e.norm()) << endl;

    // Exporting files for LIS calculations
    saveMarketVector(b, "./b.mtx");

    return 0;
}