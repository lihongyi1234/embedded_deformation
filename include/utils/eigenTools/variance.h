#ifndef EIGENTOOLS_VARIANCE_H
#define EIGENTOOLS_VARIANCE_H

#include <Eigen/Core>

inline double variance(Eigen::VectorXd& input) {
  int input_size = std::max(input.cols(), input.rows());
  Eigen::VectorXd input_mean =
      Eigen::VectorXd::Constant(input_size, input.mean());
  return (input - input_mean).array().square().sum() / input.cols();
}

#endif
