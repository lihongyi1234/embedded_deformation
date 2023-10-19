#ifndef GETMINMAX_H
#define GETMINMAX_H

#include <Eigen/Core>

static inline void getMinMax(Eigen::MatrixXd& in_cloud,
                             Eigen::Vector3d& min_point,
                             Eigen::Vector3d& max_point) {
  if (in_cloud.rows() == 3) {
    max_point = in_cloud.rowwise().maxCoeff();
    min_point = in_cloud.rowwise().minCoeff();
  } else if (in_cloud.cols() == 3) {
    max_point = in_cloud.rowwise().maxCoeff();
    min_point = in_cloud.rowwise().minCoeff();
  } else {
    throw std::invalid_argument("Error in " + std::string(__func__) + ": wrong input size");
  }
};

inline void getScale(Eigen::MatrixXd in_cloud, double& scale) {
  Eigen::Vector3d min_point;
  Eigen::Vector3d max_point;

  getMinMax(in_cloud, min_point, max_point);

  scale = (max_point - min_point).norm();
};

#endif