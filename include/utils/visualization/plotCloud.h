/**
 * Author: R. Falque
 *
 * plot cloud in polyscope
 * by R. Falque
 * 26/09/2019
 **/

#ifndef PLOT_CLOUD_HPP
#define PLOT_CLOUD_HPP

#include <Eigen/Core>
#include <string>

#include "polyscope/point_cloud.h"
#include "polyscope/view.h"

#ifndef POLYSCOPE_IS_INITIALIZED
#define POLYSCOPE_IS_INITIALIZED
int polyscope_is_initialized{0};
#endif

// polyscope wrapper
class CloudVisualization {
 private:
  std::string cloud_object_name_ = "cloud";

 public:
  CloudVisualization() { init(); }

  // destructor
  ~CloudVisualization() { polyscope::removeAllStructures(); }

  void init() {
    if (!polyscope_is_initialized) {
      polyscope::init();
      polyscope_is_initialized = 1;
    }

    // Options
    polyscope::options::autocenterStructures = false;
    polyscope::view::style = polyscope::view::NavigateStyle::Free;
    polyscope::options::groundPlaneMode = polyscope::GroundPlaneMode::None;
    polyscope::view::upDir = polyscope::view::UpDir::ZUp;
    polyscope::view::windowWidth = 1024;
    polyscope::view::windowHeight = 1024;
  }

  void add_cloud(const Eigen::MatrixXd& cloud) {
    polyscope::registerPointCloud(cloud_object_name_, cloud.transpose());
    polyscope::getPointCloud(cloud_object_name_)
        ->setPointColor(glm::vec3{0.1, 0.1, 1});
    polyscope::getPointCloud(cloud_object_name_)->setMaterial("flat");
    polyscope::getPointCloud(cloud_object_name_)->setPointRadius(0.00005, true);
    polyscope::view::resetCameraToHomeView();
  }

  void add_color(const Eigen::MatrixXd& color, std::string color_name) {
    if ((color.array() > 1.0).any() || (color.array() < 0.0).any())
      throw std::invalid_argument(
          "Error in " + std::string(__func__) + ": color contais elements out "
          "of bounds (it should be between 0 and 1)");
    if (color.rows() != 0) {
      polyscope::getPointCloud(cloud_object_name_)
          ->addColorQuantity(color_name, color.transpose());
      polyscope::getPointCloud(cloud_object_name_)
          ->getQuantity(color_name)
          ->setEnabled(true);
    }
  }

  void screenshot(std::string screenshot_path) {
    polyscope::screenshot(screenshot_path, false);
  }

  void show() { polyscope::show(); }
};

inline bool plot_cloud(const Eigen::MatrixXd& V) {
  CloudVisualization viz;
  viz.add_cloud(V);
  viz.show();
  return true;
};

inline bool plot_cloud(const Eigen::MatrixXd& V, const Eigen::MatrixXd& color) {
  CloudVisualization viz;
  viz.add_cloud(V);
  viz.add_color(color, "highlight");
  viz.show();
  return true;
};

inline bool plot_cloud(const Eigen::MatrixXd& V, const Eigen::MatrixXd& color,
                       std::string screenshot_path) {
  CloudVisualization viz;
  viz.add_cloud(V);
  viz.add_color(color, "highlight");
  viz.screenshot(screenshot_path);
  viz.show();
  return true;
};

inline bool screenshot_cloud(const Eigen::MatrixXd& V,
                             const Eigen::MatrixXd& color,
                             std::string screenshot_path) {
  CloudVisualization viz;
  viz.add_cloud(V);
  viz.add_color(color, "highlight");
  viz.screenshot(screenshot_path);
  return true;
};

#endif
