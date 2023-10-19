// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_LIST_TO_MATRIX_H
#define IGL_LIST_TO_MATRIX_H
#include <Eigen/Core>
#include <cassert>
#include <cstdio>
#include <vector>

#include "max_size.h"
#include "min_size.h"

namespace igl {
template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<std::vector<T> >& V,
                           Eigen::PlainObjectBase<Derived>& M);
template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<std::vector<T> >& V, const int n,
                           const T& padding,
                           Eigen::PlainObjectBase<Derived>& M);
template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<T>& V,
                           Eigen::PlainObjectBase<Derived>& M);

// Convert a list (std::vector) of row vectors of the same length to a matrix
// Template:
//   T  type that can be safely cast to type in Mat via '='
//   Mat  Matrix type, must implement:
//     .resize(m,n)
//     .row(i) = Row
// Inputs:
//   V  a m-long list of vectors of size n
// Outputs:
//   M  an m by n matrix
// Returns true on success, false on errors

template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<std::vector<T> >& V,
                           Eigen::PlainObjectBase<Derived>& M) {
  // number of rows
  int m = V.size();
  if (m == 0) {
    M.resize(0, 0);
    return true;
  }
  // number of columns
  int n = igl::min_size(V);
  if (n != igl::max_size(V)) {
    return false;
  }
  assert(n != -1);
  // Resize output
  M.resize(m, n);

  // Loop over rows
  for (int i = 0; i < m; i++) {
    // Loop over cols
    for (int j = 0; j < n; j++) {
      M(i, j) = V[i][j];
    }
  }

  return true;
}

template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<std::vector<T> >& V, const int n,
                           const T& padding,
                           Eigen::PlainObjectBase<Derived>& M) {
  const int m = V.size();
  M.resize(m, n);
  for (int i = 0; i < m; i++) {
    const auto& row = V[i];
    if (row.size() > n) {
      return false;
    }
    int j = 0;
    for (; j < row.size(); j++) {
      M(i, j) = row[j];
    }
    for (; j < n; j++) {
      M(i, j) = padding;
    }
  }
  return true;
}

// Convert a list of row vectors of `n` or less to a matrix and pad on
// the right with `padding`:
//
// Inputs:
//   V  a m-long list of vectors of size <=n
//   n  number of columns
//   padding  value to fill in from right for short rows
// Outputs:
//   M  an m by n matrix
template <typename T, typename Derived>
inline bool list_to_matrix(const std::vector<T>& V,
                           Eigen::PlainObjectBase<Derived>& M) {
  // number of rows
  int m = V.size();
  if (m == 0) {
    // fprintf(stderr,"Error: list_to_matrix() list is empty()\n");
    // return false;
    if (Derived::ColsAtCompileTime == 1) {
      M.resize(0, 1);
    } else if (Derived::RowsAtCompileTime == 1) {
      M.resize(1, 0);
    } else {
      M.resize(0, 0);
    }
    return true;
  }
  // Resize output
  if (Derived::RowsAtCompileTime == 1) {
    M.resize(1, m);
  } else {
    M.resize(m, 1);
  }

  // Loop over rows
  for (int i = 0; i < m; i++) {
    M(i) = V[i];
  }

  return true;
}

}  // namespace igl

#endif
