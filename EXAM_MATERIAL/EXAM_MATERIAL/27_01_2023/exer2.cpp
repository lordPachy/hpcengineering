#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;
using SpMat = SparseMatrix<double>;

int main(){
    /*******************************************************************************************************
     * PART 6
     *******************************************************************************************************/
    // Building the matrix A
     constexpr int n = 200;
    SpMat A(n, n);
    for (int i = 1; i < n - 1; i++){
        A.coeffRef(i, i-1) = -i+0.;
        A.coeffRef(i, i+1) = -i-1.;
        A.coeffRef(i, i) = 2.*(i+1.);
    }

    A.coeffRef(0, 0) = 2.;
    A.coeffRef(0, 1) = -1.;
    A.coeffRef(n-1, n-2) = -(n-1.);
    A.coeffRef(n-1, n-1) = 2.*(n);

    SpMat B(n, n);
    for (int i = 0; i < n; i++){
        B.coeffRef(i, n - 1 - i) = 1;
    }

    A += B;

    cout << "The norm of A is " << A.norm() << endl;

    /*******************************************************************************************************
     * PART 7
     *******************************************************************************************************/
    MatrixXd A_dense = MatrixXd(A);
    SelfAdjointEigenSolver<MatrixXd> eigensolver(A_dense);
    if (eigensolver.info() != Eigen::Success) abort();
    std::cout << "The first two eigenvalues of A are:\n"
        << eigensolver.eigenvalues()[0] << "\n"
        << eigensolver.eigenvalues()[1] 
        << std::endl;
 
    /*******************************************************************************************************
     * PART 8
     *******************************************************************************************************/
    saveMarket(A, "exer2.mtx");
       
    return 0;
}