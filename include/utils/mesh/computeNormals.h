/*
 *   compute the normals using the faces, can handle rowise or colwise inputs
 *   by R. Falque
 *   16/01/2020
 */

#ifndef COMPUTE_NORMALS_H
#define COMPUTE_NORMALS_H

#include <Eigen/Dense>
#include <iostream>

Eigen::MatrixXd compute_vertices_normals(Eigen::MatrixXd V, Eigen::MatrixXi F);
Eigen::MatrixXd compute_faces_normals(Eigen::MatrixXd V, Eigen::MatrixXi F);

#endif
