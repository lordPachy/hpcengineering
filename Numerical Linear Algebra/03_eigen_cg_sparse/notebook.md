
## 3. The CG Sparse iterative solver of Eigen

Eigen provides a built-in `Eigen::ConjugateGradient` solver. This class allows to solve for $A\boldsymbol{x} = \boldsymbol{b}$ linear problems using an iterative conjugate gradient algorithm. The matrix A must be selfadjoint. The matrix A and the vectors x and b can be either dense or sparse.

This class follows the sparse solver concept and has the following inputs:
- `MatrixType_`	the type of the matrix $A$, can be a dense or a sparse matrix.
- `UpLo_` the triangular part that will be used for the computations. It can be Lower, Upper (especially used in the symmetric case), or Lower|Upper in which the full matrix entries will be considered. 
- `Preconditioner_` the type of the preconditioner. Default is DiagonalPreconditioner which has $$D = diag(A)$$ and $$D^{-1}Ax = D^{-1}b$$ which will slightly better up the performance.

The maximal number of iterations and tolerance value can be controlled via the setMaxIterations() and setTolerance() methods. The defaults are the size of the problem for the maximal number of iterations and NumTraits<Scalar>::epsilon() for the tolerance.

The tolerance corresponds to the relative residual error: 
$$
tol = |A\boldsymbol{x}- \boldsymbol{b}|/|\boldsymbol{b}|
$$

N.B. Even though the default value of `UpLo_` is `Lower`, significantly higher performance is achieved when using a complete matrix and `Lower|Upper` as the `UpLo_` template parameter.

### 3.1 Example

- Download the matrix `Asym.mtx` from webeep Lab2 folder and move it to the working directory.
- Display the size of the matrix and check if it is symmetric. 
- Take as exact solution a vector `xe` defined as in the previous example and compute the right-hand side `b`. 
- Solve the resulting linear system using the Conjugate Gradient (CG) solver available in Eigen. 
- Compute and display the relative error between the exact solution `xe` and the approximated solution.

```
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
    Eigen::loadMarket(mat, "Asym.mtx");
    
    // Check matrix properties
    std::cout << "Matrix size:"<< mat.rows() << "X" << mat.cols() << endl;
    std::cout << "Non zero entries:" << mat.nonZeros() << endl;
    SpMat B = SpMat(mat.transpose()) - mat;  // Check symmetry
    std::cout << "Norm of skew-symmetric part: " << B.norm() << endl;

    // Create Rhs b
    SpVec e = SpVec::Ones(mat.rows());    // Define exact solution
    SpVec b = mat*e;                      // Compute rhs
    SpVec x(mat.rows());

    // Set parameters for solver
    double tol = 1.e-8;                 // Convergence tolerance
    int result, maxit = 1000;           // Maximum iterations
    Eigen::DiagonalPreconditioner<double> D(mat); // Create diag preconditioner

    // Solving 
    Eigen::ConjugateGradient<SpMat, Eigen::Lower|Eigen::Upper> cg;
    cg.setMaxIterations(maxit);
    cg.setTolerance(tol);
    cg.compute(mat);
    x = cg.solve(b);
    std::cout << " Eigen native CG" << endl;
    std::cout << "#iterations:     " << cg.iterations() << endl;
    std::cout << "relative residual: " << cg.error()      << endl;
    std::cout << "effective error: " << (x-e).norm() << endl;

    return 0;    
}
```

### 3.2 Exercise

- In Eigen, construct the $50\times 50$ symmetric matrix $A$ defined such that
$$ 
A = \begin{pmatrix}
    2 & -1 & 0 & 0&\ldots & 0  \\
    -1 & 2 & -1 & 0& \ldots & 0  \\
    0 & -1 & \ddots  & \ddots &\ldots  & \vdots \\
    0 & 0 & \ddots  & \ddots  & \ddots & 0 \\
   \vdots& \vdots &  \vdots &\ddots &\ddots  & -1\\
    0 & 0  &\ldots & 0& -1   & 2
\end{pmatrix}.
$$
- Display the number of nonzero entries and check if it is symmetric. 
- Take as exact solution a vector `xe` defined as in the previous example and compute the right-hand side `b`. 
- Solve the resulting linear system using the Conjugate Gradient (CG) solver available in Eigen. 
- Compute and display the relative error between the exact solution `xe` and the approximated solution.
- Export matrix $A$ and the right-hand side vector $\boldsymbol{b}$ in the `.mtx` format. Move the files in the test folder of the LIS library and repeat the previous exercise using LIS.

### HOMEWORK
Repeat the previous exercise using the BiCGSTAB solver of Eigen.