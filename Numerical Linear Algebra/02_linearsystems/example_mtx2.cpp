#include <iostream>
#include <Eigen/Sparse>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

/**
 * In this example, the error should well be under machine precision (10^(-18)).
 * Note that denormal (lower than machine precision numbers) can be represented, however
 * it is not possible to treat them for calculation. The definition is
 * 1 + eps = 1
 */

int main(int argc, char** argv){
    // Load matrix
    SparseMatrix<double> mat;
    loadMarket(mat, "mhd416a.mtx");
    VectorXd b = VectorXd::Zero(mat.rows());
    loadMarketVector(b, "b.mtx");

    VectorXd xe = VectorXd::Constant(mat.rows(), 1);      
    VectorXd b2 = mat*xe; // compute new right-hand side
    cout << "difference = " << (b-b2).norm() << endl;
    return 0;    
}
