/*
 *   fill holes to form a closed surface - wrapper towards the CGAL
 * implementation by R. Falque 06/04/2021
 */

#ifndef FILL_HOLES_H
#define FILL_HOLES_H

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>

#include <Eigen/Dense>
#include <iostream>

#include "utils/convertors/igl_to_cgal.h"

inline void fill_holes(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  // convert to CGAL format
  CGAL_Poly CGAL_mesh;
  
  convert_mesh_from_IGL_to_CGAL_Polygon(V, F, CGAL_mesh);

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for (CGAL_Poly_Halfedge_handle h : halfedges(CGAL_mesh)) {
    if (h->is_border()) {
      std::vector<CGAL_Poly_Facet_handle> patch_facets;
      std::vector<CGAL_Poly_Vertex_handle> patch_vertices;
      // bool success = std::get<0>(
      // CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
      CGAL::Polygon_mesh_processing::triangulate_and_refine_hole(
          CGAL_mesh, h, std::back_inserter(patch_facets),
          std::back_inserter(patch_vertices),
          CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
              get(CGAL::vertex_point, CGAL_mesh))
              .geom_traits(Kernel()));
      // std::cout << " Number of facets in constructed patch: " <<
      // patch_facets.size() << std::endl; std::cout << " Number of vertices in
      // constructed patch: " << patch_vertices.size() << std::endl; std::cout
      // << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
      ++nb_holes;
    }
  }

  // get back into eigen format
  convert_mesh_from_CGAL_poly_to_IGL(CGAL_mesh, V, F);
}

inline void fill_holes_and_fair(Eigen::MatrixXd& V, Eigen::MatrixXi& F) {
  // convert to CGAL format
  CGAL_Poly CGAL_mesh;
  
  convert_mesh_from_IGL_to_CGAL_Polygon(V, F, CGAL_mesh);

  // Incrementally fill the holes
  unsigned int nb_holes = 0;
  for (CGAL_Poly_Halfedge_handle h : halfedges(CGAL_mesh)) {
    if (h->is_border()) {
      std::vector<CGAL_Poly_Facet_handle> patch_facets;
      std::vector<CGAL_Poly_Vertex_handle> patch_vertices;
      // bool success = std::get<0>(
      // CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
      CGAL::Polygon_mesh_processing::triangulate_refine_and_fair_hole(
          CGAL_mesh, h, std::back_inserter(patch_facets),
          std::back_inserter(patch_vertices),
          CGAL::Polygon_mesh_processing::parameters::vertex_point_map(
              get(CGAL::vertex_point, CGAL_mesh))
              .geom_traits(Kernel()));
      // std::cout << " Number of facets in constructed patch: " <<
      // patch_facets.size() << std::endl; std::cout << " Number of vertices in
      // constructed patch: " << patch_vertices.size() << std::endl; std::cout
      // << " Fairing : " << (success ? "succeeded" : "failed") << std::endl;
      ++nb_holes;
    }
  }

  // get back into eigen format
  convert_mesh_from_CGAL_poly_to_IGL(CGAL_mesh, V, F);
}

#endif
