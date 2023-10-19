// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_WRITEOFF_H
#define IGL_WRITEOFF_H

#include <Eigen/Core>
#include <cstdio>
#include <fstream>
#include <string>

namespace igl {
// Export geometry and colors-by-vertex
// Export a mesh from an ascii OFF file, filling in vertex positions.
// Only triangle meshes are supported
//
// Templates:
//   Scalar  type for positions and vectors (will be read as double and cast
//     to Scalar)
//   Index  type for indices (will be read as int and cast to Index)
// Inputs:
//  str  path to .off output file
//   V  #V by 3 mesh vertex positions
//   F  #F by 3 mesh indices into V
//   C  double matrix of rgb values per vertex #V by 3
// Outputs:
// Returns true on success, false on errors
template <typename DerivedV, typename DerivedF, typename DerivedC>
inline bool writeOFF(const std::string str,
                     const Eigen::PlainObjectBase<DerivedV>& V,
                     const Eigen::PlainObjectBase<DerivedF>& F,
                     const Eigen::PlainObjectBase<DerivedC>& C);

template <typename DerivedV, typename DerivedF>
inline bool writeOFF(const std::string str,
                     const Eigen::PlainObjectBase<DerivedV>& V,
                     const Eigen::PlainObjectBase<DerivedF>& F);
}  // namespace igl

// write mesh to an ascii off file
template <typename DerivedV, typename DerivedF>
inline bool igl::writeOFF(const std::string fname,
                          const Eigen::PlainObjectBase<DerivedV>& V,
                          const Eigen::PlainObjectBase<DerivedF>& F) {
  using namespace std;
  using namespace Eigen;
  assert(V.cols() == 3 && "V should have 3 columns");
  ofstream s(fname);
  if (!s.is_open()) {
    fprintf(stderr, "IOError: writeOFF() could not open %s\n", fname.c_str());
    return false;
  }

  s << "OFF\n"
    << V.rows() << " " << F.rows() << " 0\n"
    << V.format(
           IOFormat(FullPrecision, DontAlignCols, " ", "\n", "", "", "", "\n"))
    << (F.array()).format(IOFormat(FullPrecision, DontAlignCols, " ", "\n",
                                   "3 ", "", "", "\n"));
  return true;
}

// write mesh and colors-by-vertex to an ascii off file
template <typename DerivedV, typename DerivedF, typename DerivedC>
inline bool igl::writeOFF(const std::string fname,
                          const Eigen::PlainObjectBase<DerivedV>& V,
                          const Eigen::PlainObjectBase<DerivedF>& F,
                          const Eigen::PlainObjectBase<DerivedC>& C) {
  using namespace std;
  using namespace Eigen;
  assert(V.cols() == 3 && "V should have 3 columns");
  assert(C.cols() == 3 && "C should have 3 columns");

  if (V.rows() != C.rows()) {
    fprintf(stderr,
            "IOError: writeOFF() Only color per vertex supported. V and C "
            "should have same size.\n");
    return false;
  }

  ofstream s(fname);
  if (!s.is_open()) {
    fprintf(stderr, "IOError: writeOFF() could not open %s\n", fname.c_str());
    return false;
  }

  // Check if RGB values are in the range [0..1] or [0..255]
  int rgbScale = (C.maxCoeff() <= 1.0) ? 255 : 1;
  // Use RGB_Array instead of RGB because of clash with mingw macro
  // (https://github.com/libigl/libigl/pull/679)
  Eigen::MatrixXd RGB_Array = rgbScale * C;

  s << "COFF\n" << V.rows() << " " << F.rows() << " 0\n";
  for (unsigned i = 0; i < V.rows(); i++) {
    s << V.row(i).format(
        IOFormat(FullPrecision, DontAlignCols, " ", " ", "", "", "", " "));
    s << unsigned(RGB_Array(i, 0)) << " " << unsigned(RGB_Array(i, 1)) << " "
      << unsigned(RGB_Array(i, 2)) << " 255\n";
  }

  s << (F.array()).format(
      IOFormat(FullPrecision, DontAlignCols, " ", "\n", "3 ", "", "", "\n"));
  return true;
}

#endif
