#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>


using namespace std;
using namespace Eigen;
using SpMat = SparseMatrix<double>;
using SpVec = VectorXd;

int main(){
    SpMat A;
    SpMat B;
    SpMat C;
    loadMarket(A, "A_block.mtx");
    loadMarket(B, "B_block.mtx");
    loadMarket(C, "C_block_sparse.mtx");

    MatrixXd M_dense(A.rows() + B.rows(), A.cols() + B.rows());
    M_dense.topLeftCorner(A.rows(), A.cols()) = A;
    M_dense.topRightCorner(B.cols(), B.rows()) = B.transpose();
    M_dense.bottomLeftCorner(B.rows(), B.cols()) = B;
    M_dense.bottomRightCorner(C.rows(), C.cols()) = C;

    SpMat M(M_dense.sparseView());

    cout << "Non-zero entries of M: " << M.nonZeros() << endl;
    cout << "Total entries of M: " << M.rows() * M.cols() << endl;

    // Decomposing the matrix M
    SpVec x(M.cols());
    SpVec b = SpVec::Ones(M.rows());
    Eigen::SparseLU<SpMat> solvelu;        // LU factorization
    solvelu.compute(M);
    if(solvelu.info()!=Eigen::Success) {                          // sanity check
        std::cout << "cannot factorize the matrix" << std::endl;  // decomposition failed
        return 0;
    }

    x = solvelu.solve(b);                                         // solving
    std::cout << "Solution with Eigen LU:" << std::endl;
    std::cout << "Residual of the solution: " << (b - M*x).norm() << endl;

    // Computing Schur's complement
    Eigen::SparseLU<SpMat> solvelu_a;        // LU factorization
    solvelu_a.compute(A);
    if(solvelu_a.info()!=Eigen::Success) {                          // sanity check
        std::cout << "cannot factorize the matrix" << std::endl;  // decomposition failed
        return 0;
    }

    SpMat Id = MatrixXd::Identity(A.rows(), A.cols()).sparseView();
    SpMat A_inv = solvelu_a.solve(Id);
    SpMat B_t = B.transpose();
    SpMat S = C - B * A_inv * B_t;

    cout << "Size of S: " << S.rows() << "x" << S.cols() << endl;
    cout << "Norm of S: " << S.norm() << endl;

    return 0;
}