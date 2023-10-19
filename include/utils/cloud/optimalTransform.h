#ifndef OPTIMAL_TRANSFORM_H
#define OPTIMAL_TRANSFORM_H

#include <Eigen/Dense>
#include <Eigen/SVD>

// following http://nghiaho.com/?page_id=671
inline void optimal_transform(Eigen::MatrixXd A, Eigen::MatrixXd B,
                              Eigen::Matrix3d &R, Eigen::Vector3d &t) {
  // check input size
  if (A.rows() != 3 || B.rows() != 3 || A.cols() != B.cols() || A.cols()<3)
    throw std::invalid_argument(
        "Error in " + std::string(__func__) + ": wrong input size");

  // get centroids
  Eigen::Vector3d centroid_A = A.rowwise().mean();
  Eigen::Vector3d centroid_B = B.rowwise().mean();
  
  // SVD
  Eigen::MatrixXd H = ( A.colwise() - centroid_A ) * ( B.colwise() - centroid_B ).transpose();
  Eigen::JacobiSVD<Eigen::MatrixXd> svd(H, Eigen::ComputeThinU | Eigen::ComputeThinV);
  Eigen::Matrix3d U = svd.matrixU();
  Eigen::Matrix3d V = svd.matrixV();

  // get R
  R = V * U.transpose();
  if (R.determinant() < 0) {
      V.col(2) = -V.col(2);
      R = V * U.transpose();
  }

  // get t
  t = centroid_B - R*centroid_A;
};

#endif
