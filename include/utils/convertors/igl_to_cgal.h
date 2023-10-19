#ifndef IGL_TO_CGAL
#define IGL_TO_CGAL

/*
 *   convertion between IGL and CGAL
 */

#include <Eigen/Core>

#include "cgal_typedef.h"

void convert_mesh_from_IGL_to_CGAL(const Eigen::MatrixXd& V,
                                   const Eigen::MatrixXi& F,
                                   CGAL_Mesh& CGAL_mesh);
void convert_mesh_from_IGL_to_CGAL_Polygon(const Eigen::MatrixXd& V,
                                           const Eigen::MatrixXi& F,
                                           CGAL_Poly& CGAL_mesh);

void convert_mesh_from_CGAL_to_IGL(const CGAL_Mesh& CGAL_mesh,
                                   Eigen::MatrixXd& V, Eigen::MatrixXi& F);
void convert_mesh_from_CGAL_poly_to_IGL(const CGAL_Poly& CGAL_mesh,
                                        Eigen::MatrixXd& V, Eigen::MatrixXi& F);

#endif
