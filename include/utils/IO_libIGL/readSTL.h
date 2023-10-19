// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READSTL_H
#define IGL_READSTL_H

#ifndef IGL_NO_EIGEN
#include <Eigen/Core>
#endif
#include <cstdio>
#include <iostream>
#include <string>
#include <vector>

#include "list_to_matrix.h"

namespace igl {
template <typename DerivedV, typename DerivedF, typename DerivedN>
inline bool readSTL(const std::string& filename,
                    Eigen::PlainObjectBase<DerivedV>& V,
                    Eigen::PlainObjectBase<DerivedF>& F,
                    Eigen::PlainObjectBase<DerivedN>& N);
template <typename TypeV, typename TypeF, typename TypeN>
inline bool readSTL(const std::string& filename,
                    std::vector<std::vector<TypeV> >& V,
                    std::vector<std::vector<TypeF> >& F,
                    std::vector<std::vector<TypeN> >& N);
template <typename TypeV, typename TypeF, typename TypeN>
inline bool readSTL(FILE* stl_file, std::vector<std::vector<TypeV> >& V,
                    std::vector<std::vector<TypeF> >& F,
                    std::vector<std::vector<TypeN> >& N);

// Read a mesh from an ascii/binary stl file.
//
// Templates:
//   Scalar  type for positions and vectors (will be read as double and cast
//     to Scalar)
// Inputs:
//   filename path to .stl file
// Outputs:
//   V  double matrix of vertex positions  #V*3 by 3
//   F  index matrix of triangle indices #F by 3
//   N  double matrix of surface normals #F by 3
// Returns true on success, false on errors
//
// Example:
//   bool success = readSTL(filename,temp_V,F,N);
//   remove_duplicate_vertices(temp_V,0,V,SVI,SVJ);
//   for_each(F.data(),F.data()+F.size(),[&SVJ](int & f){f=SVJ(f);});
//   writeOBJ("Downloads/cat.obj",V,F);
template <typename DerivedV, typename DerivedF, typename DerivedN>
inline bool readSTL(const std::string& filename,
                    Eigen::PlainObjectBase<DerivedV>& V,
                    Eigen::PlainObjectBase<DerivedF>& F,
                    Eigen::PlainObjectBase<DerivedN>& N) {
  using namespace std;
  vector<vector<typename DerivedV::Scalar> > vV;
  vector<vector<typename DerivedN::Scalar> > vN;
  vector<vector<typename DerivedF::Scalar> > vF;
  if (!readSTL(filename, vV, vF, vN)) {
    return false;
  }

  if (!list_to_matrix(vV, V)) {
    return false;
  }

  if (!list_to_matrix(vF, F)) {
    return false;
  }

  if (!list_to_matrix(vN, N)) {
    return false;
  }
  return true;
}

template <typename TypeV, typename TypeF, typename TypeN>
inline bool readSTL(const std::string& filename,
                    std::vector<std::vector<TypeV> >& V,
                    std::vector<std::vector<TypeF> >& F,
                    std::vector<std::vector<TypeN> >& N) {
  using namespace std;
  // Should test for ascii

  // Open file, and check for error
  FILE* stl_file = fopen(filename.c_str(), "rb");
  if (NULL == stl_file) {
    fprintf(stderr, "IOError: %s could not be opened...\n", filename.c_str());
    return false;
  }
  return readSTL(stl_file, V, F, N);
}

template <typename TypeV, typename TypeF, typename TypeN>
inline bool readSTL(FILE* stl_file, std::vector<std::vector<TypeV> >& V,
                    std::vector<std::vector<TypeF> >& F,
                    std::vector<std::vector<TypeN> >& N) {
  using namespace std;
  // stl_file = freopen(NULL,"rb",stl_file);
  if (NULL == stl_file) {
    fprintf(stderr,
            "IOError: stl file could not be reopened as binary (1) ...\n");
    return false;
  }

  V.clear();
  F.clear();
  N.clear();

  // Specifically 80 character header
  char header[80];
  char solid[80];
  bool is_ascii = true;
  if (fread(header, 1, 80, stl_file) != 80) {
    cerr << "IOError: too short (1)." << endl;
    goto close_false;
  }
  sscanf(header, "%s", solid);
  if (string("solid") != solid) {
    // definitely **not** ascii
    is_ascii = false;
  } else {
    // might still be binary
    char buf[4];
    if (fread(buf, 1, 4, stl_file) != 4) {
      cerr << "IOError: too short (3)." << endl;
      goto close_false;
    }
    size_t num_faces = *reinterpret_cast<unsigned int*>(buf);
    fseek(stl_file, 0, SEEK_END);
    int file_size = ftell(stl_file);
    if (file_size == 80 + 4 + (4 * 12 + 2) * num_faces) {
      is_ascii = false;
    } else {
      is_ascii = true;
    }
  }

  if (is_ascii) {
    // Rewind to end of header
    // stl_file = fopen(filename.c_str(),"r");
    // stl_file = freopen(NULL,"r",stl_file);
    fseek(stl_file, 0, SEEK_SET);
    if (NULL == stl_file) {
      fprintf(stderr, "IOError: stl file could not be reopened as ascii ...\n");
      return false;
    }
    // Read 80 header
    // Eat file name
#ifndef IGL_LINE_MAX
#define IGL_LINE_MAX 2048
#endif
    char name[IGL_LINE_MAX];
    if (NULL == fgets(name, IGL_LINE_MAX, stl_file)) {
      cerr << "IOError: ascii too short (2)." << endl;
      goto close_false;
    }
    // ascii
    while (true) {
      int ret;
      char facet[IGL_LINE_MAX], normal[IGL_LINE_MAX];
      vector<TypeN> n(3);
      double nd[3];
      ret = fscanf(stl_file, "%s %s %lg %lg %lg", facet, normal, nd, nd + 1,
                   nd + 2);
      if (string("endsolid") == facet) {
        break;
      }
      if (ret != 5 || !(string("facet") == facet || string("faced") == facet) ||
          string("normal") != normal) {
        cout << "facet: " << facet << endl;
        cout << "normal: " << normal << endl;
        cerr << "IOError: bad format (1)." << endl;
        goto close_false;
      }
      // copy casts to Type
      n[0] = nd[0];
      n[1] = nd[1];
      n[2] = nd[2];
      N.push_back(n);
      char outer[IGL_LINE_MAX], loop[IGL_LINE_MAX];
      ret = fscanf(stl_file, "%s %s", outer, loop);
      if (ret != 2 || string("outer") != outer || string("loop") != loop) {
        cerr << "IOError: bad format (2)." << endl;
        goto close_false;
      }
      vector<TypeF> f;
      while (true) {
        char word[IGL_LINE_MAX];
        int ret = fscanf(stl_file, "%s", word);
        if (ret == 1 && string("endloop") == word) {
          break;
        } else if (ret == 1 && string("vertex") == word) {
          vector<TypeV> v(3);
          double vd[3];
          int ret = fscanf(stl_file, "%lg %lg %lg", vd, vd + 1, vd + 2);
          if (ret != 3) {
            cerr << "IOError: bad format (3)." << endl;
            goto close_false;
          }
          f.push_back(V.size());
          // copy casts to Type
          v[0] = vd[0];
          v[1] = vd[1];
          v[2] = vd[2];
          V.push_back(v);
        } else {
          cerr << "IOError: bad format (4)." << endl;
          goto close_false;
        }
      }
      F.push_back(f);
      char endfacet[IGL_LINE_MAX];
      ret = fscanf(stl_file, "%s", endfacet);
      if (ret != 1 || string("endfacet") != endfacet) {
        cerr << "IOError: bad format (5)." << endl;
        goto close_false;
      }
    }
    // read endfacet
    goto close_true;
  } else {
    // Binary
    // stl_file = freopen(NULL,"rb",stl_file);
    fseek(stl_file, 0, SEEK_SET);
    // Read 80 header
    char header[80];
    if (fread(header, sizeof(char), 80, stl_file) != 80) {
      cerr << "IOError: bad format (6)." << endl;
      goto close_false;
    }
    // Read number of triangles
    unsigned int num_tri;
    if (fread(&num_tri, sizeof(unsigned int), 1, stl_file) != 1) {
      cerr << "IOError: bad format (7)." << endl;
      goto close_false;
    }
    V.resize(num_tri * 3, vector<TypeV>(3, 0));
    N.resize(num_tri, vector<TypeN>(3, 0));
    F.resize(num_tri, vector<TypeF>(3, 0));
    for (int t = 0; t < (int)num_tri; t++) {
      // Read normal
      float n[3];
      if (fread(n, sizeof(float), 3, stl_file) != 3) {
        cerr << "IOError: bad format (8)." << endl;
        goto close_false;
      }
      // Read each vertex
      for (int c = 0; c < 3; c++) {
        F[t][c] = 3 * t + c;
        N[t][c] = n[c];
        float v[3];
        if (fread(v, sizeof(float), 3, stl_file) != 3) {
          cerr << "IOError: bad format (9)." << endl;
          goto close_false;
        }
        V[3 * t + c][0] = v[0];
        V[3 * t + c][1] = v[1];
        V[3 * t + c][2] = v[2];
      }
      // Read attribute size
      unsigned short att_count;
      if (fread(&att_count, sizeof(unsigned short), 1, stl_file) != 1) {
        cerr << "IOError: bad format (10)." << endl;
        goto close_false;
      }
    }
    goto close_true;
  }
close_false:
  fclose(stl_file);
  return false;
close_true:
  fclose(stl_file);
  return true;
}
}  // namespace igl

#endif
