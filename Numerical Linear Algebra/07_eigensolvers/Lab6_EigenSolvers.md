# Eigenvalues and Eigenvectors solvers with LIS and Eigen

Lis (Library of Iterative Solvers for linear systems) provides a set of iterative methods to solve the Eigenvalue problem $A\boldsymbol{x} = \lambda \boldsymbol{x}$. 

As a first example, we aim to run the test `etest1` presented in the user guide. To do so, we follow the next steps:

- Enable `mk` modules by tiping 
```
source /u/sw/etc/bash.bashrc
``` 

- Load the LIS module: `module load lis`

- Move into the test folder: `cd lis-2.1.6/test/`

- Compile the file `etest1.c` by typing 

```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest1.c -o eigen1
```

- Run the code with multiprocessors using mpi by typing

```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt 
```
In this case, the residual is defined as $$r := (A - \lambda I) v$$

## Options for eigensolvers 

Several additional options can be selected. A non-exhaustive list is presented in the following examples:

- to chose the eigensolver `-e + [method]` (`pi` for power method that returns the largest eigenvalue, while `ii`, `rqi`, and `cr` give the smallest eigenvalue):

```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e pi
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e cr
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e rqi
```

- for the methods resorting to inner linear iterative solvers, we can specify which solvers we want to use by the option `-i + [method]` (cg, bicg, gmres, bicgstab, gs, ...). They can also be preconditioned by the option `-p + [method]` (same methods available for sparse linear system):

```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii -i cg
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii -i gs -p ilu 
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e ii -i gmres
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e rqi -i bicgstab -p ssor
```
Note that using a preconditioner does not change the number of iterations required to perform the method; it may make the inner iterative solver faster. In the same way, the inner iterative solver may make the overall method faster; the number of iterations will be unaltered.

- As for iterative methods for sparse linear system, we can specify the maximum number of iterations and the desired tolerance by setting the options `-emaxiter + [int]` and `-etol + [real]`: 

```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e pi -emaxiter 100 -etol 1.e-6
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii -emaxiter 200 -etol 1.e-15
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e ii -i gmres -etol 1.e-14
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e ii -i gs -etol 1.e-14
```

- The option `-shift + [real]` can be used to accelerate the convergence or to compute eigenpairs different from the one corresponding to the largest or smallest eigenvalues of $A$. Given the shift $\mu$, the selected method is applied to the matrix $A -\mu I_d$, where $I_d$ is the identity matrix of the same size of $A$.
The convergence of the P.M. depends on $\frac{\lambda_n}{\lambda_{n-1}}$. A shift could make this faster, essentially. 
```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e pi -shift 4.0
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii -shift 2.0
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e ii -i bicgstab -shift 8.0
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e ii -shift 1.0
```
It is relevant to note that a $\mu$ of high module could change the absolute value of the eigenvalues in such a way that the biggest eigenvalue in module is different. This may be useful in order to discover, for instance, the smallest one; however, one must be careful about the desired output.
## Shifted inverse PM
The formula goes as $$(A-\mu I)^{-1}v = \frac{1}{\lambda - \mu}v$$ The computed eigenvalue will be the closest to $\mu$.\
It is possible that picking a $\mu=\lambda$, we will get a non-invertible matrix. However, this event has probability 0; the real problem is that getting too close to some $\lambda_i$ will make the matrix $(A-\mu I)$ close to singular, thus making the "inversion" lengthy. It is a trade-off between getting close to that, in such a way the acceleration is higher, and the inversion problem $(A-\mu I)^{-1}$ is easy.
## Compute multiple eigenpairs

In LIS, some numerical methods to compute multiple eigenpairs are also available. The choice of the particular method can be specified by the option `-e + [method]`, while the number of computed eigenvalues (also referred to as $\textit{the size of the subspace}$) by the option `-ss + [int]`. This methods hinges on computing matrices similar to $A$ in Hessenberg or tridiagonal forms by the application of Arnoldi or Lanczos procedures. 

A subspace method for eigenvalue problems basically works as follows:

- A basis for a search space $span(V)$ is expanded in each step 

- An appropriate approximate eigenpair is extracted from the search subspace in each step 

- An appropriate solution of the projected problem is lifted to an approximate solution in $n$-space. 

- Convergence is monitored via the stopping criterion.

Methods differ in the way the expansion vector is selected. In LIS the methods are named after the expansion strategy. As first examples, we run `etest1` with the following options:

```
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e si -ss 4 
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e li -ss 4 -ie cg
mpirun -n 4 ./eigen1 testmat0.mtx eigvec.txt hist.txt -e li -ss 4 -ie ii -i bicgstab -p jacobi
mpirun -n 4 ./eigen1 testmat2.mtx eigvec.txt hist.txt -e ai -ss 2 -ie rqi
```

Another implemented solution for computing multiple eigenpairs in LIS is given by `etest5`. To assess the example, follow the instructions below.

- Compile the file `etest5.c` by typing 

```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest5.c -o eigen2
```

- Run the code with multiprocessors using mpi by typing

```
mpirun -n 4 ./eigen2 testmat0.mtx evals.mtx eigvecs.mtx res.txt iters.txt -ss 4 -e li 

mpirun -n 4 ./eigen2 testmat0.mtx  evals.mtx eigvecs.mtx res.txt iters.txt -ss 4 -e li -i cg -p jacobi -etol 1.0e-10 

mpirun -n 4 ./eigen2 testmat0.mtx evals.mtx eigvecs.mtx r.txt iters.txt -e si -ie ii -ss 4 -i cg -p ssor -etol 1.0e-8

mpirun -n 4 ./eigen2 testmat2.mtx evals.mtx eigvecs.mtx res.txt iters.txt -e ai -si ii -i gmres -ss 4
```


## Exercise 1

- 1. Download the symmetric matrix [fluid_sym.mtx]. Compute the largest eigenvalue of the matrix up to a tolerance of order $10^{-6}$.

- 2.  Compute the smallest eigenvalue of `fluid_sym.mtx` using the Inverse method. Explore different iterative methods and preconditioners in order to achieve a precision smaller than $10^{-10}$. Compare and comment the results.

- 3. Compute the eigenvalue of `fluid_sym.mtx` closest to different positive value of $\mu$ using the Inverse method with shift. Explore different iterative methods and comment the results.

- 4. Simultaneously compute four eigenvalues of the matrix and save the corresponding eigenvectors in a $\texttt{.mtx}$ file. Set a precision of magnitude $10^{-7}$. 


## Exercise 2

- 1. Following the instructions available on the LIS user-guide, compile and run `etest2.c` and `etest4.c`

```
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest2.c -o etest2
mpicc -DUSE_MPI -I${mkLisInc} -L${mkLisLib} -llis etest4.c -o etest4

mpirun -n 4 ./etest2 20 20 1 eigvec.mtx hist.txt
mpirun -n 4 ./etest2 20 20 1 eigvec.mtx hist.txt -e pi
mpirun -n 4 ./etest2 20 20 1 eigvec.mtx hist.txt -e rqi -i gmres
mpirun -n 4 ./etest2 20 20 1 eigvec.mtx hist.txt -e ii -i gs
mpirun -n 4 ./etest2 20 20 1 eigvec.mtx hist.txt -e ii -i cg -p ssor

mpirun -n 4 ./etest4 100
mpirun -n 4 ./etest4 100 -e pi -etol 1.0e-8
mpirun -n 4 ./etest4 50 -e pi -etol 1.0e-8 -emaxiter 2000
mpirun -n 4 ./etest4 50 -e ii -etol 1.0e-10 -i cg
mpirun -n 4 ./etest4 50 -e ii -etol 1.0e-10 -i bicgstab -p jacobi
```

- 2. Solve the eigenproblem considered in `etest4.c` using Eigen assuming $n=20, 50, 100$.

```
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
    int n = 50;	
    SparseMatrix<double> mat(n,n);            // define matrix
    for (int i=0; i<n; i++) {
        mat.coeffRef(i, i) = 2.0;
	      if(i>0) mat.coeffRef(i, i-1) = -1.0;
        if(i<n-1) mat.coeffRef(i, i+1) = -1.0;	
    }

   MatrixXd A;
   A = MatrixXd(mat);
   SelfAdjointEigenSolver<MatrixXd> eigensolver(A);
   if (eigensolver.info() != Eigen::Success) abort();
   std::cout << "The eigenvalues of A are:\n" << eigensolver.eigenvalues() << std::endl;
   // std::cout << "Here's a matrix whose columns are eigenvectors of A \n"
   //           << eigensolver.eigenvectors() << std::endl;
    return 0;    
}
```

## Exercise 3

Download and unzip the matrix `nos1.mtx` from the [matrix market website](https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/lanpro/nos1.mtx.gz). 
```
wget https://math.nist.gov/pub/MatrixMarket2/Harwell-Boeing/lanpro/nos1.mtx.gz
gzip -dk nos1.mtx.gz
```

Load the matrix in a `.cpp` file and compute its eigenvalues using the `EigenSolver` function available in Eigen. Compute the eigenvalues of the symmetric part of the previous matrix using the `SelfAdjointEigenSolver`. This methods are based on decompositions, such as $QR$ or $U\Sigma V^*$.

```
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <unsupported/Eigen/SparseExtra>

using namespace std;
using namespace Eigen;

int main(int argc, char** argv)
{
  // Load matrix
  ...
  // Check matrix properties
  ...

  // Compute Eigenvalues of the original matrix
  // Eigensolver uses dense matrix format!!!
  MatrixXd A;
  A = MatrixXd(mat);    
  ...

  // Compute Eigenvalues of symmetric matrix
  ...
  return 0;
}
```
Since we are not dealing with symmetric matrices, it is possible to get complex eigenvalues; thus, they will be represented with pairs of values.
## Exercise 4

- 1. Download the symmetric matrix [stokes_sym.mtx](https://webeep.polimi.it/mod/folder/view.php?id=129876). 

- 2. Load the matrix in a `.cpp` file and compute its eigenvalues using the `SelfAdjointEigenSolver` function available in Eigen. 

- 3. Modify the loaded matrix by adding a tridiagonal matrix corresponding to the stiffness matrix arising from a finite difference discretization of a 1D Laplacian (namely having the components on the main diagonal equal to $2$ and the other coefficients on the other two diagonals equal to $-1$).

- 4. Solve again the eigenvalue problem $Ax = \lambda x$, with $A$ obtained as in the previous point using again the `SelfAdjointEigenSolver` function.

- 5. Export the matrix $A$ in a `.mtx` file in order to be loaded on LIS codes.

- 6. Using LIS compute the largest and smallest eigenvalues and compare this results with the one previously obtained on Eigen. 