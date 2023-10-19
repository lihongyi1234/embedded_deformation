/**
 * Author: R. Falque
 * 
 * plot mesh in polyscope
 * by R. Falque
 * 26/09/2019
 **/

#ifndef POLYSCOPE_WRAPPER_HPP
#define POLYSCOPE_WRAPPER_HPP

#include <Eigen/Core>

#include <string>

#include "utils/IO/structures.h"

#include "polyscope/polyscope.h"
#include "polyscope/messages.h"
#include "polyscope/surface_mesh.h"
#include "polyscope/point_cloud.h"
#include "polyscope/curve_network.h"
#include "polyscope/view.h"

namespace visualization {

    inline void init() {
            polyscope::init();
            polyscope::options::autocenterStructures = false;
            polyscope::view::windowWidth = 1024;
            polyscope::view::windowHeight = 1024;
            polyscope::view::style = polyscope::view::NavigateStyle::Free;
            polyscope::view::upDir = polyscope::view::UpDir::NegZUp;
    }

    inline void close() {
        polyscope::removeAllStructures();
    }

    // meshes
    inline void add_mesh(const Eigen::MatrixXd& vertices, const Eigen::MatrixXi& faces, std::string mesh_name) {
        polyscope::registerSurfaceMesh(mesh_name, vertices.transpose(), faces.transpose());
        polyscope::getSurfaceMesh(mesh_name)->setSurfaceColor(glm::vec3{0.1, 0.1, 1});
        polyscope::view::resetCameraToHomeView();
    }

    inline void add_colors_to_mesh(const Eigen::MatrixXd & colors, std::string mesh_name, std::string color_name) {
        if (colors.rows() != 0){
            polyscope::getSurfaceMesh(mesh_name)->addVertexColorQuantity(color_name, colors.transpose());
            polyscope::getSurfaceMesh(mesh_name)->getQuantity(color_name)->setEnabled(true);
        }
    }

    inline void add_colors_to_mesh(const Eigen::MatrixXi & colors, std::string mesh_name, std::string color_name) {
        Eigen::MatrixXd color_double = colors.cast <double> ();
        color_double = color_double/255;
        if (color_double.rows() != 0){
            polyscope::getSurfaceMesh(mesh_name)->addVertexColorQuantity(color_name, color_double.transpose());
            polyscope::getSurfaceMesh(mesh_name)->getQuantity(color_name)->setEnabled(true);
        }
    }

    inline void set_mesh_color(const Eigen::Vector3d& color, std::string mesh_name) {
        polyscope::getSurfaceMesh(mesh_name)->setSurfaceColor(glm::vec3{color(0), color(1), color(2)});
    }
    
    
    // cloud
    inline void add_cloud(const Eigen::MatrixXd & cloud, std::string cloud_name) {
        polyscope::registerPointCloud(cloud_name, cloud.transpose());
        polyscope::getPointCloud(cloud_name)->setPointColor(glm::vec3{0.1, 0.1, 1});
        polyscope::view::resetCameraToHomeView();
    }

    inline void add_colors_to_cloud(const Eigen::MatrixXd & colors, std::string cloud_name, std::string color_name) {
        if (colors.rows() != 0){
            polyscope::getPointCloud(cloud_name)->addColorQuantity(color_name, colors);
            polyscope::getPointCloud(cloud_name)->getQuantity(color_name)->setEnabled(true);
        }
    }


    // shape
    inline void add_shape(const Mesh & shape, std::string shape_name) {
        if (shape.F.cols() == 0) {
            polyscope::registerPointCloud(shape_name, shape.V.transpose());
            polyscope::getPointCloud(shape_name)->setPointColor(glm::vec3{0.1, 0.1, 1});
            polyscope::view::resetCameraToHomeView();
            if (shape.RGB.cols() != 0){
                Eigen::MatrixXd C = shape.RGB.cast <double> ()/255;
                polyscope::getPointCloud(shape_name)->addColorQuantity("color", C.transpose());
                polyscope::getPointCloud(shape_name)->getQuantity("color")->setEnabled(true);
            }
        } else {
            polyscope::registerSurfaceMesh(shape_name, shape.V.transpose(), shape.F.transpose());
            polyscope::getSurfaceMesh(shape_name)->setSurfaceColor(glm::vec3{0.1, 0.1, 1});
            polyscope::view::resetCameraToHomeView();
            if (shape.RGB.cols() != 0){
                Eigen::MatrixXd C = shape.RGB.cast <double> ()/255;
                polyscope::getSurfaceMesh(shape_name)->addVertexColorQuantity("color", C.transpose());
                polyscope::getSurfaceMesh(shape_name)->getQuantity("color")->setEnabled(true);
            }
        }
    }


	inline void add_color_to_graph(const Eigen::MatrixXd & nodes_colors,
							const Eigen::MatrixXd & edges_colors,
							std::string graph_name,
							std::string color_name) {
		if (nodes_colors.rows() != 0){
			polyscope::getPointCloud(graph_name+"_nodes")->addColorQuantity(color_name, nodes_colors);
			polyscope::getPointCloud(graph_name+"_nodes")->getQuantity(color_name)->setEnabled(true);
		}
		if (edges_colors.rows() != 0){
			polyscope::getCurveNetwork(graph_name+"_edges")->addEdgeColorQuantity(color_name, edges_colors);
			polyscope::getCurveNetwork(graph_name+"_edges")->getQuantity(color_name)->setEnabled(true);
		}
	}

	inline void add_vector_quantity_to_graph(const Eigen::VectorXd & nodes_colors,
									  const Eigen::VectorXd & edges_colors,
									  std::string graph_name,
									  std::string vector_quantity_name) {
		if (nodes_colors.rows() != 0){
			polyscope::getPointCloud(graph_name+"_nodes")->addVectorQuantity(vector_quantity_name, nodes_colors);
			polyscope::getPointCloud(graph_name+"_nodes")->getQuantity(vector_quantity_name)->setEnabled(true);
		}
		if (edges_colors.rows() != 0){
			polyscope::getCurveNetwork(graph_name+"_edges")->addEdgeVectorQuantity(vector_quantity_name, edges_colors);
			polyscope::getCurveNetwork(graph_name+"_edges")->getQuantity(vector_quantity_name)->setEnabled(true);
		}
	}

	inline void add_scalar_quantity_to_graph(const Eigen::VectorXd & nodes_colors,
									  const Eigen::VectorXd & edges_colors,
									  std::string graph_name,
									  std::string scalar_quantity_name) {
		if (nodes_colors.rows() != 0){
			polyscope::getPointCloud(graph_name+"_nodes")->addScalarQuantity(scalar_quantity_name, nodes_colors);
			polyscope::getPointCloud(graph_name+"_nodes")->getQuantity(scalar_quantity_name)->setEnabled(true);
		}
		if (edges_colors.rows() != 0){
			polyscope::getCurveNetwork(graph_name+"_edges")->addEdgeScalarQuantity(scalar_quantity_name, edges_colors);
			polyscope::getCurveNetwork(graph_name+"_edges")->getQuantity(scalar_quantity_name)->setEnabled(true);
		}
	}


    // vectors
    inline void add_vectors(const Eigen::MatrixXd& vectors_begin, const Eigen::MatrixXd& vectors_end, std::string vectors_name) {
        polyscope::registerPointCloud(vectors_name, vectors_begin.transpose());
        polyscope::getPointCloud(vectors_name)->addVectorQuantity("vectors", (vectors_end - vectors_begin).transpose(), polyscope::VectorType::STANDARD);
        polyscope::getPointCloud(vectors_name)->getQuantity("vectors")->setEnabled(true);
        polyscope::view::resetCameraToHomeView();
    }


    // tools
    inline void show() {
        polyscope::show();
    }

    inline void screenshot(std::string screenshot_path) {
        polyscope::screenshot(screenshot_path, false);
    }

    inline void clear_structures() {
        polyscope::removeAllStructures();
    }

};

#endif