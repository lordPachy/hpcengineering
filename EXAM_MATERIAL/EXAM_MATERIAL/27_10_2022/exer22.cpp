#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;
using SpMat = SparseMatrix<double>;

int main() {
    /**
     * PART 4
     */
    // Defining A
    int n = 99;
    SpMat A(99, 99);
    for (int i = 1; i < n - 1; i++){
        A.coeffRef(i, i) = abs((n + 1.)/2. - (i+1)) + 1;
        A.coeffRef(i, i-1) = 0.5;
        A.coeffRef(i, i+1) = 0.5;
    }

    A.coeffRef(0, 0) = abs((n + 1.)/2. - (0+1)) + 1;
    A.coeffRef(n-1, n-1) = abs((n + 1.)/2. - (n-1+1)) + 1;
    A.coeffRef(0, 1) = 0.5;
    A.coeffRef(n-1, n-2) = 0.5;

    // Printing a1 and an
    cout << "a_1 = " << A.coeffRef(0, 0) << endl;
    cout << "a_n = " << A.coeffRef(n-1, n-1) << endl;

    /**
     * PART 5
     */
    // Solving the eigenvalue problem
    MatrixXd A_dense;
    A_dense = MatrixXd(A);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A_dense);
    if (eigensolver.info() != Eigen::Success) abort();

    // Reporting the largest and smallest eigenvalues
    cout << "The smallest eigenvalue of A is: " << eigensolver.eigenvalues().coeffRef(0) << endl;
    cout << "The largest eigenvalue of A is: " << eigensolver.eigenvalues().coeffRef(n-1) << endl;

    /**
     * PART 6
     */
    // Saving the matrix so it can be elaborated with LIS
    saveMarket(A, "./exer2.mtx");

    return 0;
}
