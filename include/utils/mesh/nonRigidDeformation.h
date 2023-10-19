/*
 *   ARAP - wrapper towards the CGAL implementation
 *   by R. Falque
 *   17/11/2022
 */

#ifndef NONRIGIDDEFORMATION_H
#define NONRIGIDDEFORMATION_H

#include <vector>
#include <Eigen/Dense>

#include "utils/convertors/igl_to_cgal.h"

// built upon the sample codes from https://doc.cgal.org/latest/Surface_mesh_deformation/index.html
void deform(CGAL_Poly & CGAL_mesh, std::vector<int> handles_ids, std::vector<Eigen::Vector3d> handles_target);
void deform(Eigen::MatrixXd& V, Eigen::MatrixXi& F, std::vector<Eigen::Vector3d> handle_sources, std::vector<Eigen::Vector3d> handles_target); 

#endif
