#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;

// Some useful alias
using SpMat=Eigen::SparseMatrix<double>;
using SpVec=Eigen::VectorXd;

int main(int argc, char** argv)
{
    // Load matrix
    SpMat mat;
    Eigen::loadMarket(mat, "./data/Asym.mtx");
    
    // Check matrix properties
    std::cout << "Matrix size:"<< mat.rows() << "X" << mat.cols() << endl;
    std::cout << "Non zero entries:" << mat.nonZeros() << endl;
    SpMat B = SpMat(mat.transpose()) - mat;  // Check symmetry: it needs to be "effectively 0", not 0
    std::cout << "Norm of skew-symmetric part: " << B.norm() << endl;

    // Create Rhs b
    SpVec e = SpVec::Ones(mat.rows());    // Define exact solution
    SpVec b = mat*e;                      // Compute rhs
    SpVec x(mat.rows());

    // Set parameters for solver
    double tol = 1.e-14;                 // Convergence tolerance
    int result, maxit = 50;           // Maximum iterations
    Eigen::DiagonalPreconditioner<double> D(mat); // Create diag preconditioner

    // Solving 
    Eigen::ConjugateGradient<SpMat, Eigen::Lower> cg;
    cg.setMaxIterations(maxit);
    cg.setTolerance(tol);
    cg.compute(mat);
    x = cg.solve(b);
    std::cout << " Eigen native CG" << endl;
    std::cout << "#iterations:     " << cg.iterations() << endl;
    // THE FUNCTION .error() DOES NOT RETURN THE ERROR!!!
    std::cout << "relative residual: " << cg.error()      << endl;
    std::cout << "effective error: " << (x-e).norm() << endl;

    return 0;    
}