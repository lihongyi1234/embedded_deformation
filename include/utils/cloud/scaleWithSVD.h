#ifndef SVD_SCALING_H
#define SVD_SCALING_H

#include <iostream>
#include <vector>
#include <Eigen/Core>
#include <Eigen/SVD>

inline void scale_with_SVD(Eigen::MatrixXd& src_V_in, const Eigen::MatrixXd &trg_V_in) { 
    Eigen::MatrixXd src_V = src_V_in.transpose();
    Eigen::MatrixXd trg_V = trg_V_in.transpose();

    // perform SVD for each input
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_source(src_V, Eigen::ComputeThinU | Eigen::ComputeThinV);
    Eigen::JacobiSVD<Eigen::MatrixXd> svd_target(trg_V, Eigen::ComputeThinU | Eigen::ComputeThinV);

    // change of basis (towards the basis defined by V) https://en.wikipedia.org/wiki/Change_of_basis
    src_V = src_V * svd_source.matrixV().inverse();

    // scaling in the V basis
    //src_V = src_V.rowwise() * svd_target.singularValues() / svd_source.singularValues();
    src_V.col(0) *= svd_target.singularValues()(0) / svd_source.singularValues()(0);
    src_V.col(1) *= svd_target.singularValues()(1) / svd_source.singularValues()(1);
    src_V.col(2) *= svd_target.singularValues()(2) / svd_source.singularValues()(2);

    // change back to XYZ basis
    src_V = src_V * svd_source.matrixV();

    src_V_in = src_V.transpose();
};

#endif
