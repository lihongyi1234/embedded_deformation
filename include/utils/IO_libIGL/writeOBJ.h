// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITEOBJ_H
#define IGL_WRITEOBJ_H
// History:
//  return type changed from void to bool  Alec 20 Sept 2011

#include <Eigen/Core>
#include <cassert>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

namespace igl {
// Write a mesh in an ascii obj file
// Inputs:
//   str  path to outputfile
//   V  #V by 3 mesh vertex positions
//   F  #F by 3|4 mesh indices into V
//   CN #CN by 3 normal vectors
//   FN  #F by 3|4 corner normal indices into CN
//   TC  #TC by 2|3 texture coordinates
//   FTC #F by 3|4 corner texture coord indices into TC
// Returns true on success, false on error
//
// Known issues: Horrifyingly, this does not have the same order of
// parameters as readOBJ.
template <typename DerivedV, typename DerivedF, typename DerivedCN,
          typename DerivedFN, typename DerivedTC, typename DerivedFTC>
inline bool writeOBJ(const std::string str,
                     const Eigen::MatrixBase<DerivedV>& V,
                     const Eigen::MatrixBase<DerivedF>& F,
                     const Eigen::MatrixBase<DerivedCN>& CN,
                     const Eigen::MatrixBase<DerivedFN>& FN,
                     const Eigen::MatrixBase<DerivedTC>& TC,
                     const Eigen::MatrixBase<DerivedFTC>& FTC);
template <typename DerivedV, typename DerivedF>
inline bool writeOBJ(const std::string str,
                     const Eigen::MatrixBase<DerivedV>& V,
                     const Eigen::MatrixBase<DerivedF>& F);

// Write a mesh of mixed tris and quads to an ascii obj file
// Inputs:
//   str  path to outputfile
//   V  #V by 3 mesh vertex positions
//   F  #F std::vector of std::vector<Index> of size 3 or 4 mesh indices into V
// Returns true on success, false on error
template <typename DerivedV, typename T>
inline bool writeOBJ(const std::string& str,
                     const Eigen::MatrixBase<DerivedV>& V,
                     const std::vector<std::vector<T> >& F);

}  // namespace igl

template <typename DerivedV, typename DerivedF, typename DerivedCN,
          typename DerivedFN, typename DerivedTC, typename DerivedFTC>
inline bool igl::writeOBJ(const std::string str,
                          const Eigen::MatrixBase<DerivedV>& V,
                          const Eigen::MatrixBase<DerivedF>& F,
                          const Eigen::MatrixBase<DerivedCN>& CN,
                          const Eigen::MatrixBase<DerivedFN>& FN,
                          const Eigen::MatrixBase<DerivedTC>& TC,
                          const Eigen::MatrixBase<DerivedFTC>& FTC) {
  FILE* obj_file = fopen(str.c_str(), "w");
  if (NULL == obj_file) {
    printf("IOError: %s could not be opened for writing...", str.c_str());
    return false;
  }
  // Loop over V
  for (int i = 0; i < (int)V.rows(); i++) {
    fprintf(obj_file, "v");
    for (int j = 0; j < (int)V.cols(); ++j) {
      fprintf(obj_file, " %0.17g", V(i, j));
    }
    fprintf(obj_file, "\n");
  }
  bool write_N = CN.rows() > 0;

  if (write_N) {
    for (int i = 0; i < (int)CN.rows(); i++) {
      fprintf(obj_file, "vn %0.17g %0.17g %0.17g\n", CN(i, 0), CN(i, 1),
              CN(i, 2));
    }
    fprintf(obj_file, "\n");
  }

  bool write_texture_coords = TC.rows() > 0;

  if (write_texture_coords) {
    for (int i = 0; i < (int)TC.rows(); i++) {
      fprintf(obj_file, "vt %0.17g %0.17g\n", TC(i, 0), TC(i, 1));
    }
    fprintf(obj_file, "\n");
  }

  // loop over F
  for (int i = 0; i < (int)F.rows(); ++i) {
    fprintf(obj_file, "f");
    for (int j = 0; j < (int)F.cols(); ++j) {
      // OBJ is 1-indexed
      fprintf(obj_file, " %u", F(i, j) + 1);

      if (write_texture_coords) fprintf(obj_file, "/%u", FTC(i, j) + 1);
      if (write_N) {
        if (write_texture_coords)
          fprintf(obj_file, "/%u", FN(i, j) + 1);
        else
          fprintf(obj_file, "//%u", FN(i, j) + 1);
      }
    }
    fprintf(obj_file, "\n");
  }
  fclose(obj_file);
  return true;
}

template <typename DerivedV, typename DerivedF>
inline bool igl::writeOBJ(const std::string str,
                          const Eigen::MatrixBase<DerivedV>& V,
                          const Eigen::MatrixBase<DerivedF>& F) {
  using namespace std;
  using namespace Eigen;
  assert(V.cols() == 3 && "V should have 3 columns");
  ofstream s(str);
  if (!s.is_open()) {
    fprintf(stderr, "IOError: writeOBJ() could not open %s\n", str.c_str());
    return false;
  }
  s << V.format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "",
                         "\n"))
    << (F.array() + 1)
           .format(IOFormat(FullPrecision, DontAlignCols, " ", "\n", "f ", "",
                            "", "\n"));
  return true;
}

template <typename DerivedV, typename T>
inline bool igl::writeOBJ(const std::string& str,
                          const Eigen::MatrixBase<DerivedV>& V,
                          const std::vector<std::vector<T> >& F) {
  using namespace std;
  using namespace Eigen;
  assert(V.cols() == 3 && "V should have 3 columns");
  ofstream s(str);
  if (!s.is_open()) {
    fprintf(stderr, "IOError: writeOBJ() could not open %s\n", str.c_str());
    return false;
  }
  s << V.format(
      IOFormat(FullPrecision, DontAlignCols, " ", "\n", "v ", "", "", "\n"));

  for (const auto& face : F) {
    int face_size = face.size();
    assert(face_size != 0);

    s << (face_size == 2 ? "l" : "f");

    for (const auto& vi : face) {
      s << " " << vi;
    }
    s << "\n";
  }
  return true;
}

#endif
