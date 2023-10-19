#ifndef POISSON_SURFACE_RECONSTRUCTION_H
#define POISSON_SURFACE_RECONSTRUCTION_H

#include <Eigen/Core>

void poisson_surface_reconstruction(Eigen::MatrixXd& V,
                                    Eigen::MatrixXd& N,
                                    Eigen::MatrixXd& V_out,
                                    Eigen::MatrixXi& F_out);

#endif
