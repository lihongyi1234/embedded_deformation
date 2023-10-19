/**
 * visualize correspondences
 *
 * by R. Falque
 * 17/01/2020
 **/

#ifndef VISUALIZATION_CORRESPONDENCES_H
#define VISUALIZATION_CORRESPONDENCES_H

#include <Eigen/Core>
#include <iostream>

#include "libGraphCpp/graph.hpp"
#include "polyscope/curve_network.h"
#include "polyscope/messages.h"
#include "polyscope/point_cloud.h"
#include "polyscope/polyscope.h"
#include "polyscope/surface_mesh.h"

#ifndef POLYSCOPE_IS_INITIALIZED
#define POLYSCOPE_IS_INITIALIZED
int polyscope_is_initialized{0};
#endif

inline void createCurveNetwork(const Eigen::MatrixXd &source,
                               const Eigen::MatrixXd &target,
                               Eigen::MatrixXd &nodes, Eigen::MatrixXi &edges) {
  nodes.resize(source.rows() + target.rows(), 3);
  edges.resize(source.rows(), 2);
  for (int i = 0; i < source.rows(); i++) {
    nodes.row(i * 2) = source.row(i);
    nodes.row(i * 2 + 1) = target.row(i);
    edges.row(i) << i * 2, i * 2 + 1;
  }
};


inline void plot_correspondences(const Eigen::MatrixXd & source_V,
                                const Eigen::MatrixXi & source_F,
                                const Eigen::MatrixXd & target_V,
                                const Eigen::MatrixXd & ED_source,
                                const Eigen::MatrixXd & ED_target)
{

    if (!polyscope_is_initialized) {
        polyscope::init();
        polyscope_is_initialized = 1;
    }

    polyscope::options::autocenterStructures = false;
    polyscope::view::windowWidth = 1024;
    polyscope::view::windowHeight = 1024;

    polyscope::registerSurfaceMesh("source", source_V.transpose(), source_F.transpose());
    polyscope::getSurfaceMesh("source")->setSurfaceColor(glm::vec3{0.1, 0.1, 1});

    polyscope::registerPointCloud("target", target_V.transpose());
    polyscope::getPointCloud("target")->setPointColor(glm::vec3{1, 0.1, 0.1});

    Eigen::MatrixXd nodes;
    Eigen::MatrixXi edges;
    createCurveNetwork(ED_source.transpose(), ED_target.transpose(), nodes, edges);

    polyscope::registerCurveNetwork("correspondences", nodes, edges);
    polyscope::getCurveNetwork("correspondences")->setColor(glm::vec3{0, 0, 0});

    polyscope::show();
    polyscope::removeAllStructures();
}


inline void plot_correspondences(const Eigen::MatrixXd &source_V,
                                 const Eigen::MatrixXi &source_F,
                                 const Eigen::MatrixXd &target_V,
                                 const Eigen::MatrixXi &target_F,
                                 const Eigen::MatrixXd &ED_source,
                                 const Eigen::MatrixXd &ED_target) {
  if (!polyscope_is_initialized) {
    polyscope::init();
    polyscope_is_initialized = 1;
  }

  polyscope::options::autocenterStructures = false;
  polyscope::view::windowWidth = 1024;
  polyscope::view::windowHeight = 1024;

  polyscope::registerSurfaceMesh("source", source_V.transpose(),
                                 source_F.transpose());
  polyscope::getSurfaceMesh("source")->setSurfaceColor(glm::vec3{0.1, 0.1, 1});

  polyscope::registerSurfaceMesh("target", target_V.transpose(),
                                 target_F.transpose());
  polyscope::getSurfaceMesh("target")->setSurfaceColor(glm::vec3{1, 0.1, 0.1});

  Eigen::MatrixXd nodes;
  Eigen::MatrixXi edges;
  createCurveNetwork(ED_source.transpose(), ED_target.transpose(), nodes,
                     edges);

  polyscope::registerCurveNetwork("correspondences", nodes, edges);
  polyscope::getCurveNetwork("correspondences")->setColor(glm::vec3{0, 0, 0});

  polyscope::show();

  polyscope::removeAllStructures();
}

inline void plot_correspondences(libgraphcpp::Graph &source,
                                 libgraphcpp::Graph &target,
                                 Eigen::VectorXd &rigidity_terms) {
  Eigen::MatrixXd source_V, target_V;
  Eigen::MatrixXi source_E, target_E;

  source_V = source.get_nodes();
  target_V = target.get_nodes();
  source_E = source.get_edges();
  target_E = target.get_edges();

  Eigen::VectorXd empty;

  Eigen::MatrixXd nodes;
  Eigen::MatrixXi edges;
  createCurveNetwork(source_V, target_V, nodes, edges);

  polyscope::registerCurveNetwork("correspondences", nodes, edges);
  polyscope::getCurveNetwork("correspondences")->setColor(glm::vec3{0, 0, 0});

  // libgraphcpp::polyscope_plot_with_color(source_V, source_E, rigidity_terms,
  // empty);
}

inline void plot_correspondences(libgraphcpp::Graph &source,
                                 libgraphcpp::Graph &target,
                                 Eigen::VectorXd &rigidity_terms,
                                 double offset) {
  Eigen::MatrixXd source_V, target_V;
  Eigen::MatrixXi source_E, target_E;

  source_V = source.get_nodes();
  target_V = target.get_nodes();

  Eigen::VectorXd offset_vector =
      Eigen::VectorXd::Constant(target_V.rows(), offset);
  target_V.col(2) += offset_vector;

  source_E = source.get_edges();
  target_E = target.get_edges();

  Eigen::VectorXd empty;

  Eigen::MatrixXd correspondences_nodes;
  Eigen::MatrixXi correspondences_edges;
  createCurveNetwork(source_V, target_V, correspondences_nodes,
                     correspondences_edges);

  polyscope::removeAllStructures();
  polyscope::registerCurveNetwork("correspondences", correspondences_nodes,
                                  correspondences_edges);
  polyscope::getCurveNetwork("correspondences")->setColor(glm::vec3{1, 0, 1});

  // libgraphcpp::add_graph_with_scalar_to_plot(source_V, source_E,
  // rigidity_terms, empty);
  // libgraphcpp::add_graph_with_scalar_to_plot(target_V, target_E,
  // rigidity_terms, empty);

  polyscope::show();
  polyscope::removeAllStructures();
}

inline void plot_correspondences(libgraphcpp::Graph &source,
                                 libgraphcpp::Graph &target,
                                 Eigen::VectorXd &rigidity_terms, double offset,
                                 int correspondence_to_remove) {
  Eigen::MatrixXd source_V, target_V;
  Eigen::MatrixXi source_E, target_E;

  source_V = source.get_nodes();
  target_V = target.get_nodes();

  Eigen::VectorXd offset_vector =
      Eigen::VectorXd::Constant(target_V.rows(), offset);
  target_V.col(2) += offset_vector;

  source_E = source.get_edges();
  target_E = target.get_edges();

  Eigen::VectorXd empty;

  Eigen::MatrixXd correspondences_nodes;
  Eigen::MatrixXi correspondences_edges;
  createCurveNetwork(source_V.row(correspondence_to_remove),
                     target_V.row(correspondence_to_remove),
                     correspondences_nodes, correspondences_edges);

  std::vector<std::array<double, 3>> edgeColor(correspondences_edges.rows());
  for (size_t i = 0; i < edgeColor.size(); i++) edgeColor[i] = {{0, 0, 1}};
  edgeColor[correspondence_to_remove] = {{1, 0, 0}};

  polyscope::removeAllStructures();
  polyscope::registerCurveNetwork("correspondences", correspondences_nodes,
                                  correspondences_edges);

  polyscope::getCurveNetwork("correspondences")->setColor(glm::vec3{1, 0, 0});

  polyscope::getCurveNetwork("correspondences")
      ->addEdgeColorQuantity("to_remove", edgeColor);
  polyscope::getCurveNetwork("correspondences")
      ->getQuantity("to_remove")
      ->setEnabled(true);

  // libgraphcpp::add_graph_with_scalar_to_plot(source_V, source_E,
  // rigidity_terms, empty);
  // libgraphcpp::add_graph_with_scalar_to_plot(target_V, target_E,
  // rigidity_terms, empty);

  polyscope::show();
  polyscope::removeAllStructures();
}

#endif
