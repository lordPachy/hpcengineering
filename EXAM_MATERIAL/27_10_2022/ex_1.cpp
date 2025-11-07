#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "grad.hpp"
#include "cg.hpp"

using namespace std;
using namespace Eigen;
using namespace LinearAlgebra;
using SpMat = SparseMatrix<double>;

int main(int argc, char** argv)
{   
    /**
     * PART 6
     */
    // Load matrix
    SpMat A;
    loadMarket(A, "./diffreact.mtx");

    // Checking the shape of the matrix
    cout << "Rows: " << A.rows() << endl;
    cout << "Columns: " << A.cols() << endl;

    // Transposing and taking the norm
    SpMat A_t = A.transpose();
    cout << "difference = " << (A - A_t).norm() << endl;

    // Finding the skew-symmetric part and showing norms of A and its s.s. part
    SpMat A_ss = 0.5 * (A - A_t);
    cout << "Norm of A: " << A.norm() << endl;
    cout << "Norm of the skew-symmetric part of A: " << A_ss.norm() << endl;

    /**
     * PART 7
     */
    VectorXd x_star = VectorXd::Constant(A.cols(), 1.0);
    VectorXd b = A*x_star;
    cout << "Norm of b: " << b.norm() << endl;

    /**
     * PART 8
     */
    // Settings for the solvers
    double tol = 1.e-8;                // Convergence tolerance
    int result, maxit = 1000;             // Maximum iterations
    int maxit_grad = 100000;
    Eigen::DiagonalPreconditioner<double> D(A); // Create diagonal preconditioner


    // Solving system and reporting results
    // Gradient
    VectorXd x = VectorXd::Constant(A.cols(), 0.0);
    result = GRAD(A, x, b, D, maxit_grad, tol);        // Solve system

    std::cout <<"Gradient Method"<< endl;
    cout << "Gradient flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    std::cout << "Error norm: "<<(x-x_star).norm()<< endl;

    // Conjugate Gradient
    x *= 0.0;
    result = CG(A, x, b, D, maxit, tol);        // Solve system

    std::cout <<"Conjugate gradient"<< endl;
    cout << "CG flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;
    std::cout << "Error norm: "<<(x-x_star).norm()<< endl;

    /**
     * PART 9
     */

    return 0;    
}
