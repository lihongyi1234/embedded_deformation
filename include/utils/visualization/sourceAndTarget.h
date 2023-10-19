/**
 * visualize correspondences
 *
 * by R. Falque
 * 17/01/2020
 **/

#ifndef VISUALIZATION_SOURCE_AND_TARGET_H
#define VISUALIZATION_SOURCE_AND_TARGET_H

#include <Eigen/Core>
#include <iostream>

#include "polyscope/messages.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#ifndef POLYSCOPE_IS_INITIALIZED
#define POLYSCOPE_IS_INITIALIZED
int polyscope_is_initialized{0};
int polyscope_object_counter{0};
#endif

inline bool plot_source_and_target(const Eigen::MatrixXd& source_V,
                                   const Eigen::MatrixXi& source_F,
                                   const Eigen::MatrixXd& target_V,
                                   const Eigen::MatrixXi& target_F) {
  if (!polyscope_is_initialized) {
    polyscope::init();
    polyscope_is_initialized = 1;
  }

  polyscope::options::autocenterStructures = false;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  polyscope::registerSurfaceMesh("source", source_V, source_F);
  polyscope::getSurfaceMesh("source")->surfaceColor = glm::vec3{0.1, 0.1, 1};

  polyscope::registerSurfaceMesh("target", target_V, target_F);
  polyscope::getSurfaceMesh("target")->surfaceColor = glm::vec3{1, 0.1, 0.1};

  polyscope::show();

  polyscope::removeAllStructures();
  return true;
};

#endif
