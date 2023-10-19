#ifndef IO_READSTL_H
#define IO_READSTL_H

#include <Eigen/Core>
#include <fstream>
#include <iostream>
#include <vector>

#include "../IO_libIGL/readSTL.h"

void readSTL(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F, Eigen::MatrixXd &N, bool verbose) {
  Eigen::MatrixXd V_temp, N_temp;
  Eigen::MatrixXi F_temp;
  igl::readSTL(filepath, V_temp, F_temp, N_temp);

  V = V_temp.transpose();
  F = F_temp.transpose();
  N = N_temp.transpose();
}

void readSTL(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F, Eigen::MatrixXd &N) {
  Eigen::MatrixXd V_temp, N_temp;
  Eigen::MatrixXi F_temp;
  igl::readSTL(filepath, V_temp, F_temp, N_temp);

  V = V_temp.transpose();
  F = F_temp.transpose();
  N = N_temp.transpose();
}

void readSTL(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXi &F) {
  Eigen::MatrixXd V_temp, N_temp;
  Eigen::MatrixXi F_temp;
  igl::readSTL(filepath, V_temp, F_temp, N_temp);

  V = V_temp.transpose();
  F = F_temp.transpose();
}

void readSTL(const std::string &filepath, Eigen::MatrixXd &V,
             Eigen::MatrixXd &N) {
  Eigen::MatrixXd V_temp, N_temp;
  Eigen::MatrixXi F_temp;
  igl::readSTL(filepath, V_temp, F_temp, N_temp);

  V = V_temp.transpose();
  N = N_temp.transpose();
}

#endif