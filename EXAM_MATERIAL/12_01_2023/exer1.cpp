#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "gmres.hpp"

using namespace std;
using namespace Eigen;
using SpMat = SparseMatrix<double>;

int main(){
    /**
     * PART 6
     */
    // Loading the matrix
    SpMat A;
    loadMarket(A, "./A_stokes.mtx");

    // Checking symmetry
    SpMat A_t;
    A_t = A.transpose();
    cout << (A - A_t).norm() << endl;

    // Checking norm of A
    cout << "Norm of A: " << A.norm() << endl;

    /**
     * PART 7
     */
    VectorXd b = VectorXd::Constant(A.cols(), 1.0);

    // Options for the solver
    double tol = 1.e-8;                // Convergence tolerance
    int result, maxit = 1000;           // Maximum iterations
    int restart = 50;                  // Restart for gmres
    DiagonalPreconditioner<double> D(A);
    VectorXd x = VectorXd::Constant(A.cols(), 0.0);

    // Solving
    result = LinearAlgebra::GMRES(A, x, b, D, maxit, maxit, tol);
    cout << "GMRES   flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    //cout << "Error:                " << (x-b).norm()<< endl;
    cout << "x_1 = " << x.coeffRef(0) << endl;

    /**
     * PART 8
     */
    // Building As
    SpMat trid(A.rows(), A.cols());
    for (int i = 1; i < A.rows() - 1; i++){
        trid.coeffRef(i, i-1) = 1.;
        trid.coeffRef(i, i) = 2.;
        trid.coeffRef(i, i+1) = -1.;
    }

    trid.coeffRef(0, 0) = 2.;
    trid.coeffRef(0, 1) = -1.;
    trid.coeffRef(A.rows() - 1, A.rows() - 2) = 1.;
    trid.coeffRef(A.rows() - 1, A.rows() - 1) = 2.;

    SpMat As = A + A.norm() * trid;

    // Computing the norm of As_s
    SpMat As_t = As.transpose();
    cout << "Norm of the symmetric part of As " << 0.5*(As + As_t).norm() << endl;

    return 0;
}