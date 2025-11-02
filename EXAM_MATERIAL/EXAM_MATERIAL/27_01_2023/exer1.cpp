#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>
#include "bcgstab.hpp"

using namespace std;
using namespace Eigen;
using SpMat = SparseMatrix<double>;

int main(){
    /*******************************************************************************************************
     * PART 9
     *******************************************************************************************************/
    // Loading matrices
    SpMat A_1, A_2;
    loadMarket(A_1, "A_1.mtx");
    loadMarket(A_2, "A_2.mtx");

    // Calculating A and its information
    SpMat A = A_1 * A_2;
    cout << "Size of A: " << A.rows() << "x" << A.cols() << endl;
    cout << "Non-zero entries of A: " << A.nonZeros() << endl;
    cout << "Norm of A: " << A.norm() << endl;

    /*******************************************************************************************************
     * PART 10
     *******************************************************************************************************/ 
    // Solver options
    double tol = 1.e-10;                // Convergence tolerance
    int result, maxit = 1000;           // Maximum iterations
    DiagonalPreconditioner<double> D(A);
    VectorXd x_slu(A.cols()), x_bicg(A.cols());
    VectorXd b = VectorXd::Constant(A.rows(), 1.);

    // SparseLU solution
    SparseLU<SpMat> sparselu;
    sparselu.compute(A);
    x_slu = sparselu.solve(b);

    // BiCGStab solution
    result = LinearAlgebra::BiCGSTAB(A, x_bicg, b, D, maxit, tol);
    cout << "BiCGSTAB   flag = " << result << endl;
    cout << "iterations performed: " << maxit << endl;
    cout << "tolerance achieved  : " << tol << endl;

    //  Informations about solution
    cout << "Norm of the SparseLU solution: " << x_slu.norm() << endl;
    cout << "Norm of the BiCGSTAB solution: " << x_bicg.norm() << endl;
    cout << "Norm of the solutions' difference: " << (x_slu - x_bicg).norm() << endl;
    

    return 0;
}