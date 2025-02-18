// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2013 Alec Jacobson <alecjacobson@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_READOFF_H
#define IGL_READOFF_H
// History:
//  return type changed from void to bool  Alec 18 Sept 2011

#include <Eigen/Core>
#include <cstdio>
#include <string>
#include <vector>

#include "list_to_matrix.h"

namespace igl {

// Read a mesh from an ascii OFF file, filling in vertex positions, normals
// and texture coordinates. Mesh may have faces of any number of degree
//
// Templates:
//   Scalar  type for positions and vectors (will be read as double and cast
//     to Scalar)
//   Index  type for indices (will be read as int and cast to Index)
// Inputs:
//  str  path to .obj file
// Outputs:
//   V  double matrix of vertex positions  #V by 3
//   F  #F list of face indices into vertex positions
//   N  list of vertex normals #V by 3
//   C  list of rgb color values per vertex #V by 3
// Returns true on success, false on errors
template <typename Scalar, typename Index>
inline bool readOFF(const std::string off_file_name,
                    std::vector<std::vector<Scalar> >& V,
                    std::vector<std::vector<Index> >& F,
                    std::vector<std::vector<Scalar> >& N,
                    std::vector<std::vector<Scalar> >& C);
// Inputs:
//   off_file  pointer to already opened .off file
// Outputs:
//   off_file  closed file
template <typename Scalar, typename Index>
inline bool readOFF(FILE* off_file, std::vector<std::vector<Scalar> >& V,
                    std::vector<std::vector<Index> >& F,
                    std::vector<std::vector<Scalar> >& N,
                    std::vector<std::vector<Scalar> >& C);

#ifndef IGL_NO_EIGEN
// read mesh from a ascii off file
// Inputs:
//   str  path to .off file
// Outputs:
//   V  eigen double matrix #V by 3
//   F  eigen int matrix #F by 3
template <typename DerivedV, typename DerivedF>
inline bool readOFF(const std::string str, Eigen::PlainObjectBase<DerivedV>& V,
                    Eigen::PlainObjectBase<DerivedF>& F);

template <typename DerivedV, typename DerivedF>
inline bool readOFF(const std::string str, Eigen::PlainObjectBase<DerivedV>& V,
                    Eigen::PlainObjectBase<DerivedF>& F,
                    Eigen::PlainObjectBase<DerivedV>& N,
                    Eigen::PlainObjectBase<DerivedV>& C);
#endif

}  // namespace igl

template <typename Scalar, typename Index>
inline bool igl::readOFF(const std::string off_file_name,
                         std::vector<std::vector<Scalar> >& V,
                         std::vector<std::vector<Index> >& F,
                         std::vector<std::vector<Scalar> >& N,
                         std::vector<std::vector<Scalar> >& C) {
  using namespace std;
  FILE* off_file = fopen(off_file_name.c_str(), "r");
  if (NULL == off_file) {
    printf("IOError: %s could not be opened...\n", off_file_name.c_str());
    return false;
  }
  return readOFF(off_file, V, F, N, C);
}

template <typename Scalar, typename Index>
inline bool igl::readOFF(FILE* off_file, std::vector<std::vector<Scalar> >& V,
                         std::vector<std::vector<Index> >& F,
                         std::vector<std::vector<Scalar> >& N,
                         std::vector<std::vector<Scalar> >& C) {
  using namespace std;
  V.clear();
  F.clear();
  N.clear();
  C.clear();

  // First line is always OFF
  char header[1000];
  const std::string OFF("OFF");
  const std::string NOFF("NOFF");
  const std::string COFF("COFF");
  if (fscanf(off_file, "%s\n", header) != 1 ||
      !(string(header).compare(0, OFF.length(), OFF) == 0 ||
        string(header).compare(0, COFF.length(), COFF) == 0 ||
        string(header).compare(0, NOFF.length(), NOFF) == 0)) {
    printf(
        "Error: readOFF() first line should be OFF or NOFF or COFF, not %s...",
        header);
    fclose(off_file);
    return false;
  }
  bool has_normals = string(header).compare(0, NOFF.length(), NOFF) == 0;
  bool has_vertexColors = string(header).compare(0, COFF.length(), COFF) == 0;
  // Second line is #vertices #faces #edges
  int number_of_vertices;
  int number_of_faces;
  int number_of_edges;
  char tic_tac_toe;
  char line[1000];
  bool still_comments = true;
  while (still_comments) {
    char* temp = fgets(line, 1000, off_file);
    still_comments = (line[0] == '#' || line[0] == '\n');
  }
  sscanf(line, "%d %d %d", &number_of_vertices, &number_of_faces,
         &number_of_edges);
  V.resize(number_of_vertices);
  if (has_normals) N.resize(number_of_vertices);
  if (has_vertexColors) C.resize(number_of_vertices);
  F.resize(number_of_faces);
  // printf("%s %d %d %d\n",(has_normals ? "NOFF" :
  // "OFF"),number_of_vertices,number_of_faces,number_of_edges);
  // Read vertices
  for (int i = 0; i < number_of_vertices;) {
    char* temp = fgets(line, 1000, off_file);
    double x, y, z, nx, ny, nz;
    if (sscanf(line, "%lg %lg %lg %lg %lg %lg", &x, &y, &z, &nx, &ny, &nz) >=
        3) {
      std::vector<Scalar> vertex;
      vertex.resize(3);
      vertex[0] = x;
      vertex[1] = y;
      vertex[2] = z;
      V[i] = vertex;

      if (has_normals) {
        std::vector<Scalar> normal;
        normal.resize(3);
        normal[0] = nx;
        normal[1] = ny;
        normal[2] = nz;
        N[i] = normal;
      }

      if (has_vertexColors) {
        C[i].resize(3);
        C[i][0] = nx / 255.0;
        C[i][1] = ny / 255.0;
        C[i][2] = nz / 255.0;
      }
      i++;
    } else if (fscanf(off_file, "%[#]", &tic_tac_toe) == 1) {
      char comment[1000];
      int temp = fscanf(off_file, "%[^\n]", comment);
    } else {
      printf("Error: bad line (%d)\n", i);
      if (feof(off_file)) {
        fclose(off_file);
        return false;
      }
    }
  }
  // Read faces
  for (int i = 0; i < number_of_faces;) {
    std::vector<Index> face;
    int valence;
    if (fscanf(off_file, "%d", &valence) == 1) {
      face.resize(valence);
      for (int j = 0; j < valence; j++) {
        int index;
        if (j < valence - 1) {
          int temp = fscanf(off_file, "%d", &index);
        } else {
          int temp = fscanf(off_file, "%d%*[^\n]", &index);
        }

        face[j] = index;
      }
      F[i] = face;
      i++;
    } else if (fscanf(off_file, "%[#]", &tic_tac_toe) == 1) {
      char comment[1000];
      int temp = fscanf(off_file, "%[^\n]", comment);
    } else {
      printf("Error: bad line\n");
      fclose(off_file);
      return false;
    }
  }
  fclose(off_file);
  return true;
}

template <typename DerivedV, typename DerivedF>
inline bool igl::readOFF(const std::string str,
                         Eigen::PlainObjectBase<DerivedV>& V,
                         Eigen::PlainObjectBase<DerivedF>& F) {
  std::vector<std::vector<double> > vV;
  std::vector<std::vector<double> > vN;
  std::vector<std::vector<int> > vF;
  std::vector<std::vector<double> > vC;
  bool success = igl::readOFF(str, vV, vF, vN, vC);
  if (!success) {
    // readOFF(str,vV,vF,vN,vC) should have already printed an error
    // message to stderr
    return false;
  }
  bool V_rect = list_to_matrix(vV, V);
  if (!V_rect) {
    // list_to_matrix(vV,V) already printed error message to std err
    return false;
  }
  bool F_rect = list_to_matrix(vF, F);
  if (!F_rect) {
    // list_to_matrix(vF,F) already printed error message to std err
    return false;
  }
  return true;
}

template <typename DerivedV, typename DerivedF>
inline bool igl::readOFF(const std::string str,
                         Eigen::PlainObjectBase<DerivedV>& V,
                         Eigen::PlainObjectBase<DerivedF>& F,
                         Eigen::PlainObjectBase<DerivedV>& N,
                         Eigen::PlainObjectBase<DerivedV>& C) {
  std::vector<std::vector<double> > vV;
  std::vector<std::vector<double> > vN;
  std::vector<std::vector<int> > vF;
  std::vector<std::vector<double> > vC;
  bool success = igl::readOFF(str, vV, vF, vN, vC);
  if (!success) {
    // readOFF(str,vV,vF,vC) should have already printed an error
    // message to stderr
    return false;
  }
  bool V_rect = list_to_matrix(vV, V);
  if (!V_rect) {
    // list_to_matrix(vV,V) already printed error message to std err
    return false;
  }
  bool F_rect = list_to_matrix(vF, F);
  if (!F_rect) {
    // list_to_matrix(vF,F) already printed error message to std err
    return false;
  }

  if (vN.size()) {
    bool N_rect = list_to_matrix(vN, N);
    if (!N_rect) {
      // list_to_matrix(vN,N) already printed error message to std err
      return false;
    }
  }

  // Warning: RGB colors will be returned in the N matrix
  if (vC.size()) {
    bool C_rect = list_to_matrix(vC, C);
    if (!C_rect) {
      // list_to_matrix(vC,N) already printed error message to std err
      return false;
    }
  }

  return true;
}

#endif
