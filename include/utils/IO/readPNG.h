// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2016 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_PNG_READ_PNG_H
#define IGL_PNG_READ_PNG_H
#include <Eigen/Core>
#include <string>
#include "stb_image.h"
#include <numeric>
#include <algorithm>


struct Image {
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> R;
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> G;
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> B;
    Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic> A;
};

// Read an image from a .png file into 4 memory buffers
//
// Input:
//  png_file  path to .png file
// Output:
//  R,G,B,A texture channels
// Returns true on success, false on failure
//
inline bool readPNG(const std::string png_file,
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& R,
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& G,
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& B,
Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>& A
) {
    int cols,rows,n;
    unsigned char *data = stbi_load(png_file.c_str(), &cols, &rows, &n, 4);
    if(data == NULL) {
        return false;
    }

    R.resize(cols,rows);
    G.resize(cols,rows);
    B.resize(cols,rows);
    A.resize(cols,rows);

    for (unsigned i=0; i<rows; ++i) {
        for (unsigned j=0; j<cols; ++j) {
            R(j,rows-1-i) = data[4*(j + cols * i) + 0];
            G(j,rows-1-i) = data[4*(j + cols * i) + 1];
            B(j,rows-1-i) = data[4*(j + cols * i) + 2];
            A(j,rows-1-i) = data[4*(j + cols * i) + 3];
        }
    }

    stbi_image_free(data);
    
    return true;
};

inline bool readPNG(const std::string png_file, Image & img) {
    return readPNG(png_file, img.R, img.G, img.B, img.A);
};

inline Image readPNG(const std::string png_file) {
    Image img;
    readPNG(png_file, img.R, img.G, img.B, img.A);
    return img;
};

// load an png file and produce a mesh with colors sampled from the RGB image input
// the dimension of the mesh is proportional to the image resolution
inline bool readPNG(const std::string png_file, Eigen::MatrixXd & V, Eigen::MatrixXi & F, Eigen::MatrixXd & RGB) {
    Image img = readPNG(png_file);
    
    int cols = img.R.cols();
    int rows = img.R.rows();
    
    // set the RGB values
    RGB.resize(3, rows*cols);
    Eigen::MatrixXd r = img.R.cast <double> (); 
    Eigen::MatrixXd g = img.G.cast <double> (); 
    Eigen::MatrixXd b = img.B.cast <double> (); 

    RGB.row(0) =  Eigen::Map<const Eigen::VectorXd>(r.data(), r.size())/255;
    RGB.row(1) =  Eigen::Map<const Eigen::VectorXd>(g.data(), g.size())/255;
    RGB.row(2) =  Eigen::Map<const Eigen::VectorXd>(b.data(), b.size())/255;

    // set the XYZ coordinates
    std::vector<int> ivec(cols * rows);
    std::iota(ivec.begin(), ivec.end(), 0);

    std::vector<int> X = ivec;
    std::for_each(X.begin(), X.end(), [&cols](int &x){ x = x % cols; });
    std::vector<double> Xd(X.begin(), X.end());

    std::vector<int> Y = ivec;
    std::for_each(Y.begin(), Y.end(), [&cols](int &y){ y = y / cols; });
    std::vector<double> Yd(Y.begin(), Y.end());

    V.resize(3, rows*cols);
    V.row(0) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Xd.data(), Xd.size());
    V.row(1) = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(Yd.data(), Yd.size());
    V.row(2).setZero();

    // set the faces
    F.resize(3, (rows-1)*(cols-1)*2);
    int counter = 0;
    for (int v=0; v<rows-1; v++) {
        for (int u=0; u<cols-1; u++) {
            F.col(counter) << u + cols*v, u+1 + cols*v, u+cols + cols*v;
            counter ++;
            F.col(counter) << u+1 + cols*v, u+cols+1 + cols*v, u+cols + cols*v;
            counter ++;
        }
    }

    return false;
}

#endif
