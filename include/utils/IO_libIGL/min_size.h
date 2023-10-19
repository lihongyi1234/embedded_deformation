// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_MIN_SIZE_H
#define IGL_MIN_SIZE_H
#include <vector>

namespace igl {
// Determine min size of lists in a vector
// Template:
//   T  some list type object that implements .size()
// Inputs:
//   V  vector of list types T
// Returns min .size() found in V, returns -1 if V is empty
template <typename T>
inline int min_size(const std::vector<T>& V) {
  int min_size = -1;
  for (typename std::vector<T>::const_iterator iter = V.begin();
       iter != V.end(); iter++) {
    int size = (int)iter->size();
    // have to handle base case
    if (min_size == -1) {
      min_size = size;
    } else {
      min_size = (min_size < size ? min_size : size);
    }
  }
  return min_size;
}
}  // namespace igl

#endif
