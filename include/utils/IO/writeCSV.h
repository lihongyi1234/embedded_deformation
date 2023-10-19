#ifndef EIGEN_TOOLS_WRITECSV_H
#define EIGEN_TOOLS_WRITECSV_H

#include <Eigen/Core>
#include <fstream>
#include <string>
#include <vector>

// using namespace Eigen;

template <typename M>
void writeCSV(std::string name, M matrix) {
  const static Eigen::IOFormat CSVFormat(Eigen::StreamPrecision,
                                         Eigen::DontAlignCols, ", ", "\n");
  std::ofstream file(name.c_str());
  file << matrix.format(CSVFormat);
  file.close();
}

template <typename M>
void writeCSV(std::string filename, M matrix, std::vector<std::string> rows_ids,
              std::vector<std::string> cols_ids) {
  // test the size of the matrix first
  // if (matrix.rows() != x_ids.size() || matrix.cols() != y_ids.size()) {
  //
  //}

  // open file
  std::ofstream myfile;
  myfile.open(filename);

  // write the x_ids
  myfile << ",";
  for (int i = 0; i < cols_ids.size(); i++) myfile << cols_ids[i] << ",";
  myfile << "\n";

  // write the content
  for (int i = 0; i < rows_ids.size(); i++) {
    myfile << rows_ids[i] << ",";
    for (int j = 0; j < cols_ids.size(); j++) myfile << matrix(i, j) << ",";
    myfile << "\n";
  }

  // close file
  myfile.close();
}

#endif