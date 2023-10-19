/*
 *   compute the normals using the faces
 *   by R. Falque
 *   16/01/2020
 */

#ifndef COMPUTE_GEOMETRICAL_PROPERTIES_H
#define COMPUTE_GEOMETRICAL_PROPERTIES_H

#include <Eigen/Dense>
#include <iostream>

double compute_volume(Eigen::MatrixXd V, Eigen::MatrixXi F, Eigen::MatrixXd N);
double compute_volume(Eigen::MatrixXd V, Eigen::MatrixXi F);
double compute_surface_area(Eigen::MatrixXd V, Eigen::MatrixXi F);

#endif
