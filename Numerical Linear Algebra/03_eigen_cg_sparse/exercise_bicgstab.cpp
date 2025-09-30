#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;

// Some useful alias
using SpMat=Eigen::SparseMatrix<double>;
using SpVec=Eigen::VectorXd;
typedef Eigen::Triplet<double> T;

int main(int argc, char** argv)
{
    // Creating matrix
    SpMat mat(50, 50);
    
    std::vector<T> tripletList;
    tripletList.reserve(148);
    for(int i=1; i<50; i++) {
      tripletList.push_back(T(i, i, 2.0));
      tripletList.push_back(T(i-1, i, -1.0));
      tripletList.push_back(T(i, i-1, -1.0));
    }

    tripletList.push_back(T(0, 0, 2.0));
    mat.setFromTriplets(tripletList.begin(), tripletList.end());
    
    
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
    double tol = 1.e-8;                 // Convergence tolerance
    int result, maxit = 1000;           // Maximum iterations
    Eigen::DiagonalPreconditioner<double> D(mat); // Create diag preconditioner

    // Solving: using ConjugateGradient with only one of the two is not faster
    // but the sparse matrix can be sparser
    Eigen::BiCGSTAB<SpMat> cg;
    cg.setMaxIterations(maxit);
    cg.setTolerance(tol);
    cg.compute(mat);
    x = cg.solve(b);
    std::cout << " Eigen native CG" << endl;
    std::cout << "#iterations:     " << cg.iterations() << endl;
    // THE FUNCTION .error() DOES NOT RETURN THE ERROR!!!
    std::cout << "relative residual: " << cg.error()      << endl;
    std::cout << "effective error: " << (x-e).norm() << endl;

    // Export vector in .mtx format
    int n = b.size();
    Eigen::saveMarket(mat, "./data/A_bicgstab.mtx");

    return 0;    
}