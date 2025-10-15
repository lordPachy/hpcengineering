#include <Eigen/Core>
#include <Eigen/SVD>
#include <iostream>
#include <string>
#include <unsupported/Eigen/SparseExtra>

int main(int argc, char *argv[]){
  using namespace Eigen;

  if(argc != 2)
    {
      std::cerr << " Usage: provide matrix filename" << std::endl;
      return 1;
    }
  std::string matrixFile(argv[1]);

  // Some useful alias
  using SpMat = SparseMatrix<double>;

  // Read matrix
  SpMat A;
  loadMarket(A, matrixFile);
  MatrixXd Ad = MatrixXd(A); // convert to dense for SVD

  // compute Eigen SVD
  Eigen::BDCSVD<Eigen::MatrixXd> svd (Ad, Eigen::ComputeThinU | Eigen::ComputeThinV);
  VectorXd W = svd.singularValues();
  std::cout << "singular values: " << W << "\n";
  MatrixXd Asvd = svd.matrixU() * svd.singularValues().asDiagonal() * svd.matrixV().transpose();
  MatrixXd diff = Asvd - Ad;
  std::cout << "diff:\n" << diff.norm() << "\n";
  return 1;
}