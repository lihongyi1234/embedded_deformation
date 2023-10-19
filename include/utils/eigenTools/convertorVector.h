#ifndef EIGENTOOLS_CONVERTOR_VECTOR_HPP
#define EIGENTOOLS_CONVERTOR_VECTOR_HPP

#include <Eigen/Core>
#include <stdexcept>
#include <vector>

template <typename T>
Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> to_eigen(
    std::vector<std::vector<T>> input) {
  Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> output;

  for (auto array : input)
    if (array.size() != input[0].size())
      throw std::invalid_argument(
          "Error in " + std::string(__func__) + ": all the sub arrays do not have the same size");

  output.resize(input.size(), input[0].size());
  for (size_t i = 0; i < input.size(); i++) {
    for (size_t j = 0; j < input[i].size(); j++) {
      output(i, j) = input[i][j];
    }
  }

  return output;
};

#endif
