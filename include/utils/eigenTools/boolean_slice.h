#ifndef BOOLEAN_SLICE_H
#define BOOLEAN_SLICE_H

#include <Eigen/Core>
#include <iostream>
#include <vector>

// Equivalent of "A = A(B == b)" from Matlab
template <typename Type_A, typename Type_B>
inline Eigen::Matrix<Type_A, Eigen::Dynamic, Eigen::Dynamic>
boolean_slice_equal(Eigen::Matrix<Type_A, Eigen::Dynamic, Eigen::Dynamic>& A,
                    Eigen::Matrix<Type_B, Eigen::Dynamic, Eigen::Dynamic>& B,
                    Eigen::Matrix<Type_B, 3, 1>& b) {
  // find B==b
  std::vector<int> index_to_keep;
  for (int i = 0; i < B.cols(); i++) {
    if (B.col(i) == b) {
      index_to_keep.push_back(i);
    }
  }

  // create temp var
  Eigen::Matrix<Type_A, Eigen::Dynamic, Eigen::Dynamic> A_temp;
  A_temp.resize(3, index_to_keep.size());

  for (int i = 0; i < index_to_keep.size(); i++) {
    A_temp.col(i) = A.col(index_to_keep.at(i));
  }

  return A_temp;
}

#endif