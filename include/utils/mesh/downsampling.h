/*
 *   downsampling - wrapper towards the CGAL implementation
 *   by R. Falque
 *   06/04/2021
 */

#ifndef DOWNSAMPLING_H
#define DOWNSAMPLING_H

#include <CGAL/Surface_mesh_simplification/Policies/Edge_collapse/Count_stop_predicate.h>
#include <CGAL/Surface_mesh_simplification/edge_collapse.h>

#include <Eigen/Dense>
#include <iostream>

#include "utils/convertors/igl_to_cgal.h"

inline void downsampling(Eigen::MatrixXd& V, Eigen::MatrixXi& F,
                         int stop_predicate) {
  // convert to CGAL format
  CGAL_Poly CGAL_mesh;
  
  convert_mesh_from_IGL_to_CGAL_Polygon(V, F, CGAL_mesh);

  // perform downsampling
  CGAL::Surface_mesh_simplification::Count_stop_predicate<CGAL_Mesh> stop(
      stop_predicate);
  int r = CGAL::Surface_mesh_simplification::edge_collapse(
      CGAL_mesh, stop,
      CGAL::parameters::vertex_index_map(
          get(CGAL::vertex_external_index, CGAL_mesh))
          .halfedge_index_map(get(CGAL::halfedge_external_index, CGAL_mesh)));

  // get back into eigen format
  convert_mesh_from_CGAL_poly_to_IGL(CGAL_mesh, V, F);
}

#endif
