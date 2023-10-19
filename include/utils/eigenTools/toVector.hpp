#ifndef EIGENTOOLS_TO_VECTOR_H
#define EIGENTOOLS_TO_VECTOR_H

#include <Eigen/Dense>
#include <vector>

using std::vector;
using Eigen::Matrix;

template <typename T, int R, int C, int N>
void matrixToVectorColWise(const Matrix<T, R, C> & matrix, vector<Matrix<T, N, 1>>& vec) {
    vec.clear();
    vec.reserve(matrix.cols());
    for (int c = 0; c < matrix.cols(); c++) {
        vec.emplace_back(matrix.col(c));
    }
}

template <typename T, int R, int C, int N>
void matrixToVectorRowWise(const Matrix<T, R, C> & matrix, vector<Matrix<T, N, 1>>& vec) {
    matrixToVectorColWise<T, C, R, N>(matrix.transpose(), vec);
}


#endif