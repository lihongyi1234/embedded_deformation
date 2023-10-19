/*
 *   Matlab style boolean tests
 *   by R. Falque
 *   12/02/2020
 */

#ifndef BOOLEAN_TESTS_H
#define BOOLEAN_TESTS_H

#include <Eigen/Dense>
#include <iostream>
#include <vector>

template <typename Type_B>
inline std::vector<int> is_equal(
    Eigen::Matrix<Type_B, Eigen::Dynamic, Eigen::Dynamic>& B,
    Eigen::Matrix<Type_B, 3, 1>& b) {
  std::vector<int> indexes;
  for (int i = 0; i < B.cols(); i++) {
    if (B.col(i) == b) {
      indexes.push_back(i);
    }
  }

  return indexes;
};

template <typename Type_B>
inline std::vector<int> is_not_equal(
    Eigen::Matrix<Type_B, Eigen::Dynamic, Eigen::Dynamic>& B,
    Eigen::Matrix<Type_B, 3, 1>& b) {
  std::vector<int> indexes;
  for (int i = 0; i < B.cols(); i++) {
    if (B.col(i) != b) {
      indexes.push_back(i);
    }
  }

  return indexes;
};

#endif
