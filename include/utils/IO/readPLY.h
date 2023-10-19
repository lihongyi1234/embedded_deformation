#ifndef IO_READPLY_H
#define IO_READPLY_H

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <vector>

#define TINYPLY_IMPLEMENTATION
#include "structures.h"
#include "tinyply.h"

using namespace tinyply;

void readPLY(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F, Eigen::MatrixXd &N, Eigen::MatrixXi &RGB,
             bool verbose) {
  try {
    std::ifstream ss(filepath, std::ios::binary);
    if (ss.fail()) throw std::runtime_error("failed to open " + filepath);

    PlyFile file;
    file.parse_header(ss);

    if (verbose) {
      std::cout << "..........................................................."
                   ".............\n";
      for (auto c : file.get_comments())
        std::cout << "Comment: " << c << std::endl;
      for (auto e : file.get_elements()) {
        std::cout << "element - " << e.name << " (" << e.size << ")"
                  << std::endl;
        for (auto p : e.properties)
          std::cout << "\tproperty - " << p.name << " ("
                    << tinyply::PropertyTable[p.propertyType].str << ")"
                    << std::endl;
      }
      std::cout << "..........................................................."
                   ".............\n";
    }
    // Tinyply treats parsed data as untyped byte buffers. See below for
    // examples.
    std::shared_ptr<tinyply::PlyData> vertices_handle, normals_handle,
        faces_handle, texcoords_handle, RGB_handle;
    bool vertices_are_float = false;
    bool normals_are_float = false;

    // get type
    for (auto e : file.get_elements()) {
      for (auto p : e.properties) {
        if (e.name == "vertex" && p.name == "x")
          vertices_are_float =
              tinyply::PropertyTable[p.propertyType].str == "float";
        if (e.name == "vertex" && p.name == "nx")
          normals_are_float =
              tinyply::PropertyTable[p.propertyType].str == "float";
      }
    }
    if (verbose) {
      std::cout << "vertices_are_float : " << vertices_are_float << std::endl;
    }

    // The header information can be used to programmatically extract properties
    // on elements known to exist in the header prior to reading the data. For
    // brevity of this sample, properties like vertex position are hard-coded:
    try {
      vertices_handle =
          file.request_properties_from_element("vertex", {"x", "y", "z"});
    } catch (const std::exception &e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    try {
      normals_handle =
          file.request_properties_from_element("vertex", {"nx", "ny", "nz"});
    } catch (const std::exception &e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    try {
      texcoords_handle =
          file.request_properties_from_element("vertex", {"u", "v"});
    } catch (const std::exception &e) {
      if (verbose) {
        std::cerr << "tinyply exception: " << e.what() << std::endl;
      }
    }

    try {
      RGB_handle = file.request_properties_from_element(
          "vertex", {"red", "green", "blue"});
    } catch (const std::exception &e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    // Providing a list size hint (the last argument) is a 2x performance
    // improvement. If you have arbitrary ply files, it is best to leave this 0.
    try {
      faces_handle =
          file.request_properties_from_element("face", {"vertex_indices"}, 3);
    } catch (const std::exception &e) {
      std::cerr << "tinyply exception: " << e.what() << std::endl;
    }

    file.read(ss);

    if (verbose) {
      if (vertices_handle)
        std::cout << "\tRead " << vertices_handle->count << " total vertices "
                  << std::endl;
      if (normals_handle)
        std::cout << "\tRead " << normals_handle->count
                  << " total vertex normals " << std::endl;
      if (faces_handle)
        std::cout << "\tRead " << faces_handle->count
                  << " total faces (triangles) " << std::endl;
      if (RGB_handle)
        std::cout << "\tRead " << RGB_handle->count << " total vertex RGB "
                  << std::endl;
    }

    if (vertices_handle) {
      if (vertices_are_float) {
        const size_t numVerticesBytes = vertices_handle->buffer.size_bytes();
        std::vector<Eigen::Vector3f> vertices(vertices_handle->count);
        std::memcpy(vertices.data(), vertices_handle->buffer.get(),
                    numVerticesBytes);
        V.resize(3, vertices.size());
        for (int i = 0; i < vertices.size(); i++)
          V.col(i) = vertices[i].cast<double>();
      } else {
        const size_t numVerticesBytes = vertices_handle->buffer.size_bytes();
        std::vector<Eigen::Vector3d> vertices(vertices_handle->count);
        std::memcpy(vertices.data(), vertices_handle->buffer.get(),
                    numVerticesBytes);
        V.resize(3, vertices.size());
        for (int i = 0; i < vertices.size(); i++) V.col(i) = vertices[i];
      }
    }

    if (normals_handle) {
      if (normals_are_float) {
        const size_t numNormalsBytes = normals_handle->buffer.size_bytes();
        std::vector<Eigen::Vector3f> normals(normals_handle->count);
        std::memcpy(normals.data(), normals_handle->buffer.get(),
                    numNormalsBytes);

        N.resize(3, normals.size());
        for (int i = 0; i < normals.size(); i++)
          N.col(i) = normals[i].cast<double>();
      } else {
        const size_t numNormalsBytes = normals_handle->buffer.size_bytes();
        std::vector<Eigen::Vector3d> normals(normals_handle->count);
        std::memcpy(normals.data(), normals_handle->buffer.get(),
                    numNormalsBytes);

        N.resize(3, normals.size());
        for (int i = 0; i < normals.size(); i++) N.col(i) = normals[i];
      }
    }

    if (faces_handle) {
      const size_t numFacesBytes = faces_handle->buffer.size_bytes();
      std::vector<Eigen::Vector3i> faces(faces_handle->count);
      std::memcpy(faces.data(), faces_handle->buffer.get(), numFacesBytes);

      F.resize(3, faces.size());
      for (int i = 0; i < faces.size(); i++) F.col(i) = faces[i];
    }
    if (RGB_handle) {
      const size_t numRgbBytes = RGB_handle->buffer.size_bytes();
      std::vector<uchar3> rgb(RGB_handle->count);
      std::memcpy(rgb.data(), RGB_handle->buffer.get(), numRgbBytes);

      RGB.resize(3, rgb.size());
      for (int i = 0; i < rgb.size(); i++)
        RGB.col(i) << int(rgb[i].r), int(rgb[i].g), int(rgb[i].b);
    }

  } catch (const std::exception &e) {
    std::cerr << "Caught tinyply exception: " << e.what() << std::endl;
  }
}

void readPLY(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F, Eigen::MatrixXd &N, Eigen::MatrixXi &RGB) {
  bool verbose = false;
  readPLY(filepath, V, F, N, RGB, verbose);
}

void readPLY(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F, Eigen::MatrixXd &N) {
  Eigen::MatrixXi RGB;
  bool verbose = false;
  readPLY(filepath, V, F, N, RGB, verbose);
}

void readPLY(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F) {
  Eigen::MatrixXd N;
  Eigen::MatrixXi RGB;
  bool verbose = false;
  readPLY(filepath, V, F, N, RGB, verbose);
}

void readPLY(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXd &N) {
  Eigen::MatrixXi F;
  Eigen::MatrixXi RGB;
  bool verbose = false;
  readPLY(filepath, V, F, N, RGB, verbose);
}


Mesh readPLY(const std::string & filepath)
{
  Mesh mesh;
  readPLY(filepath, mesh.V, mesh.F, mesh.N, mesh.RGB, false);
  return mesh;
}

#endif