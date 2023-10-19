#ifndef UNIQUE_H
#define UNIQUE_H

#include <Eigen/Core>
#include <iostream>
#include <vector>

void unique(Eigen::MatrixXi& M, std::vector<Eigen::Vector3i>& unique_elements);
void unique(Eigen::MatrixXi& M, std::vector<Eigen::Vector3i>& unique_elements,
            bool verbose);

#endif
