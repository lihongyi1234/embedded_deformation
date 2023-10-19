#ifndef EIGENTOOLS_REDUCE_H
#define EIGENTOOLS_REDUCE_H

#include <Eigen/Dense>

using Eigen::Matrix;

template <typename T, int R, int C>
T meanPairwiseDistanceColWise(const Matrix<T, R, C>& matrix) {
    T sum = 0;
    for (int i = 0; i < matrix.cols(); i++) {
        for (int j = i + 1; j < matrix.cols(); j++) {
            sum += (matrix.col(i) - matrix.col(j)).norm();
        }
    }
    return sum / static_cast<T>(matrix.cols() * (matrix.cols() - 1) / 2);
}

template <typename T, int R, int C>
T meanPairwiseDistanceRowWise(const Matrix<T, R, C>& matrix) {
    return meanPairwiseDistanceColWise<T, C, R>(matrix.transpose());
}

#endif