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
    SpMat M1(321, 321);
    SpMat M2(321, 321);

    // Main diagonal
    for (int i = 0; i < 321; i++){
        M1.coeffRef(i, i) = 4;
        M2.coeffRef(i, i) = 5;
    }

    // 1-off diagonals
    for (int i = 0+1; i < 321-1; i++){
        M1.coeffRef(i, i-1) = -3;
        M1.coeffRef(i, i+1) = -3;
        M2.coeffRef(i, i+1) = -2;
        M2.coeffRef(i, i-1) = -2;
    }
    M1.coeffRef(0, 1) = -3;
    M2.coeffRef(0, 1) = -2;
    M1.coeffRef(321-1, 321-2) = -3;
    M2.coeffRef(321-1, 321-2) = -2;

    // 2-off diagonals
    for (int i = 0+2; i < 321-2; i++){
        M1.coeffRef(i, i-2) = 2;
        M1.coeffRef(i, i+2) = 2;
    }
    M1.coeffRef(0, 2) = 2;
    M1.coeffRef(1, 3) = 2;
    M1.coeffRef(321-1, 321-3) = 2;
    M1.coeffRef(321-2, 321-4) = 2;

    // 3-off diagonals
    for (int i = 0+3; i < 321-3; i++){
        M1.coeffRef(i, i-3) = -1;
        M1.coeffRef(i, i+3) = -1;
    }
    M1.coeffRef(0, 3) = -1;
    M1.coeffRef(1, 4) = -1;
    M1.coeffRef(2, 5) = -1;
    M1.coeffRef(321-1, 321-4) = -1;
    M1.coeffRef(321-2, 321-5) = -1;
    M1.coeffRef(321-3, 321-6) = -1;

    SpMat M = M1*M2;
    SpMat M_T = M.transpose();

    cout << "Norm of skew-symmetric part of M: " << (0.5*(M-M_T)).norm() << endl;
    cout << "Norm of M: " << (M).norm() << endl;

    return 0;
}