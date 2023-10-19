#ifndef EIGENTOOLS_STANDARD_DEVIATION_H
#define EIGENTOOLS_STANDARD_DEVIATION_H

#include <Eigen/Core>

inline double standard_deviation(Eigen::VectorXd& input) {
  int input_size = std::max(input.cols(), input.rows());
  Eigen::VectorXd input_mean =
      Eigen::VectorXd::Constant(input_size, input.mean());
  return sqrt((input - input_mean).array().square().sum() / input_size);
}

#endif
