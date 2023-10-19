#pragma once

#include <ceres/ceres.h>
#include <Eigen/Dense>

// equation (3) in KinEtre: norm2(R'*R-I)^2
// gradient computed with http://www.matrixcalculus.org/

// equation (3) in KinEtre: norm2(R'*R-I)^2
struct RotCostFunction_2_V2
{
	// constructor
	RotCostFunction_2_V2(const double cost_function_weight) {
		cost_function_weight_ = sqrtf(cost_function_weight);
	}

	bool operator()(
		double const* const parameters,
		double* residuals) const
	{
		Eigen::Map<const Eigen::Matrix3d> R(parameters);
		Eigen::MatrixXd I = Eigen::MatrixXd::Identity(3, 3);

		Eigen::Map<Eigen::MatrixXd>res(residuals, 3, 3);
		res = cost_function_weight_ * (R.transpose() * R - I);
		return true;
	}

	double cost_function_weight_;
};

// equation (7)
struct RegCostFunction_V2
{
	// constructor
	RegCostFunction_V2(const double cost_function_weight, Eigen::Vector3d g_j, Eigen::Vector3d g_k) {
		g_j_ = g_j;
		g_k_ = g_k;
		cost_function_weight_ = sqrtf(cost_function_weight);
	}

	bool operator()(
		double const* const parameters_j,
		double const* const parameters_k,
		double* residuals)const
	{

		Eigen::Map<const Eigen::Matrix3d> R_j(parameters_j);
		Eigen::Map<const Eigen::Vector3d> t_j(parameters_j + 9);
		Eigen::Map<const Eigen::Matrix3d> R_k(parameters_k);
		Eigen::Map<const Eigen::Vector3d> t_k(parameters_k + 9);

		Eigen::Map<Eigen::Vector3d>res(residuals, 3);

		double alpha_j_k = 1;

		res = cost_function_weight_ * alpha_j_k * (R_j * (g_k_ - g_j_) + g_j_ + t_j - (g_k_ + t_k));

		return true;
	}

	Eigen::Vector3d g_j_;
	Eigen::Vector3d g_k_;
	double cost_function_weight_;

};

// equation (7)
struct RegCostFunction_V3
{
	// constructor
	RegCostFunction_V3(const double cost_function_weight, Eigen::Vector3d g_j, Eigen::Vector3d g_k) {
		g_j_ = g_j;
		g_k_ = g_k;
		cost_function_weight_ = sqrtf(cost_function_weight);
	}

	bool operator()(
		double const* const parameters_j,
		double const* const parameters_k,
		double* residuals)const
	{
		Eigen::Map<const Eigen::Vector3d> wxyz_j(parameters_j);
		Eigen::Quaterniond q_j(wxyz_j[0], wxyz_j[1], wxyz_j[2], wxyz_j[3]);
		q_j = q_j.normalized();
		Eigen::Matrix3d R_j = q_j.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_j(parameters_j + 4);
		Eigen::Map<const Eigen::Vector3d> t_k(parameters_k + 4);

		Eigen::Map<Eigen::Vector3d>res(residuals, 3);

		double alpha_j_k = 1;

		res = cost_function_weight_ * alpha_j_k * (R_j * (g_k_ - g_j_) + g_j_ + t_j - (g_k_ + t_k));

		return true;
	}

	Eigen::Vector3d g_j_;
	Eigen::Vector3d g_k_;
	double cost_function_weight_;

};

// equation (8)
struct ConCostFunction_V2
{
	// constructor
	ConCostFunction_V2(double cost_function_weight,
		std::vector<Eigen::Vector3d> vector_g,
		Eigen::Vector3d source,
		Eigen::Vector3d target) {
		// define residual and block_size
		source_ = source;
		target_ = target;
		vector_g_ = vector_g;
		cost_function_weight_ = sqrt(cost_function_weight);

		nodes_connectivity = vector_g_.size() - 1;

		// equation (3)
		w_j_.clear();
		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_.push_back(pow(1 - (source_ - vector_g_[i]).squaredNorm() / (source_ - vector_g_.back()).squaredNorm(), 2));

		double normalization_factor = 0;
		for (int i = 0; i < nodes_connectivity; ++i)
			normalization_factor += w_j_[i];

		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_[i] /= normalization_factor;
	}

	bool operator()(
		double const* const parameters_0,
		double const* const parameters_1,
		double const* const parameters_2,
		double const* const parameters_3,
		double const* const parameters_4,
		double const* const parameters_5,
		double* residuals) const
	{
		Eigen::Map<Eigen::MatrixXd>res(residuals, 1, 3);

		// back to equation (8)
		Eigen::Vector3d new_node_position;
		new_node_position << 0, 0, 0;
		
		/*
		for (int i = 0; i < nodes_connectivity; ++i)
		{
			Eigen::Map<const Eigen::Matrix3d> R_j(parameters + 12 * i);
			Eigen::Map<const Eigen::Vector3d> t_j(parameters + 9 + 12 * i);
			new_node_position += w_j_[i] * (R_j * (source_ - vector_g_[i]) + vector_g_[i] + t_j);
		}*/

		Eigen::Map<const Eigen::Matrix3d> R_0(parameters_0);
		Eigen::Map<const Eigen::Vector3d> t_0(parameters_0 + 9);
		new_node_position += w_j_[0] * (R_0 * (source_ - vector_g_[0]) + vector_g_[0] + t_0);

		Eigen::Map<const Eigen::Matrix3d> R_1(parameters_1);
		Eigen::Map<const Eigen::Vector3d> t_1(parameters_1 + 9);
		new_node_position += w_j_[1] * (R_1 * (source_ - vector_g_[1]) + vector_g_[1] + t_1);

		Eigen::Map<const Eigen::Matrix3d> R_2(parameters_2);
		Eigen::Map<const Eigen::Vector3d> t_2(parameters_2 + 9);
		new_node_position += w_j_[2] * (R_2 * (source_ - vector_g_[2]) + vector_g_[2] + t_2);

		Eigen::Map<const Eigen::Matrix3d> R_3(parameters_3);
		Eigen::Map<const Eigen::Vector3d> t_3(parameters_3 + 9);
		new_node_position += w_j_[3] * (R_3 * (source_ - vector_g_[3]) + vector_g_[3] + t_3);

		Eigen::Map<const Eigen::Matrix3d> R_4(parameters_4);
		Eigen::Map<const Eigen::Vector3d> t_4(parameters_4 + 9);
		new_node_position += w_j_[4] * (R_4 * (source_ - vector_g_[4]) + vector_g_[4] + t_4);

		Eigen::Map<const Eigen::Matrix3d> R_5(parameters_5);
		Eigen::Map<const Eigen::Vector3d> t_5(parameters_5 + 9);
		new_node_position += w_j_[5] * (R_5 * (source_ - vector_g_[5]) + vector_g_[5] + t_5);

		res(0, 0) = cost_function_weight_ * (new_node_position - target_)(0);
		res(0, 1) = cost_function_weight_ * (new_node_position - target_)(1);
		res(0, 2) = cost_function_weight_ * (new_node_position - target_)(2);
		return true;
	}

	Eigen::Vector3d source_;
	Eigen::Vector3d target_;
	std::vector<Eigen::Vector3d> vector_g_;
	std::vector<double> w_j_;
	int nodes_connectivity;
	double cost_function_weight_;

};


// equation (8)
struct ConCostFunction_V3
{
	// constructor
	ConCostFunction_V3(double cost_function_weight,
		std::vector<Eigen::Vector3d> vector_g,
		Eigen::Vector3d source,
		Eigen::Vector3d target) {
		// define residual and block_size
		source_ = source;
		target_ = target;
		vector_g_ = vector_g;
		cost_function_weight_ = sqrt(cost_function_weight);

		nodes_connectivity = vector_g_.size() - 1;

		// equation (3)
		w_j_.clear();
		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_.push_back(pow(1 - (source_ - vector_g_[i]).squaredNorm() / (source_ - vector_g_.back()).squaredNorm(), 2));

		double normalization_factor = 0;
		for (int i = 0; i < nodes_connectivity; ++i)
			normalization_factor += w_j_[i];

		for (int i = 0; i < nodes_connectivity; ++i)
			w_j_[i] /= normalization_factor;
	}

	bool operator()(
		double const* const parameters_0,
		double const* const parameters_1,
		double const* const parameters_2,
		double const* const parameters_3,
		double const* const parameters_4,
		double const* const parameters_5,
		double* residuals) const
	{
		Eigen::Map<Eigen::MatrixXd>res(residuals, 1, 3);

		// back to equation (8)
		Eigen::Vector3d new_node_position;
		new_node_position << 0, 0, 0;

		//Eigen::Map<const Eigen::Matrix3d> R_0(parameters_0);
		Eigen::Map<const Eigen::Vector4d> wxyz_0(parameters_0);
		Eigen::Quaterniond q_0(wxyz_0[0], wxyz_0[1], wxyz_0[2], wxyz_0[3]);
		q_0 = q_0.normalized();
		Eigen::Matrix3d R_0 = q_0.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_0(parameters_0 + 4);
		new_node_position += w_j_[0] * (R_0 * (source_ - vector_g_[0]) + vector_g_[0] + t_0);

		//Eigen::Map<const Eigen::Matrix3d> R_1(parameters_1);
		Eigen::Map<const Eigen::Vector4d> wxyz_1(parameters_1);
		Eigen::Quaterniond q_1(wxyz_1[0], wxyz_1[1], wxyz_1[2], wxyz_1[3]);
		q_1 = q_1.normalized();
		Eigen::Matrix3d R_1 = q_1.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_1(parameters_1 + 4);
		new_node_position += w_j_[1] * (R_1 * (source_ - vector_g_[1]) + vector_g_[1] + t_1);

		//Eigen::Map<const Eigen::Matrix3d> R_2(parameters_2);
		Eigen::Map<const Eigen::Vector4d> wxyz_2(parameters_2);
		Eigen::Quaterniond q_2(wxyz_2[0], wxyz_2[1], wxyz_2[2], wxyz_2[3]);
		q_2 = q_2.normalized();
		Eigen::Matrix3d R_2 = q_2.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_2(parameters_2 + 4);
		new_node_position += w_j_[2] * (R_2 * (source_ - vector_g_[2]) + vector_g_[2] + t_2);

		//Eigen::Map<const Eigen::Matrix3d> R_3(parameters_3);
		Eigen::Map<const Eigen::Vector4d> wxyz_3(parameters_3);
		Eigen::Quaterniond q_3(wxyz_3[0], wxyz_3[1], wxyz_3[2], wxyz_3[3]);
		q_3 = q_3.normalized();
		Eigen::Matrix3d R_3 = q_3.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_3(parameters_3 + 4);
		new_node_position += w_j_[3] * (R_3 * (source_ - vector_g_[3]) + vector_g_[3] + t_3);

		//Eigen::Map<const Eigen::Matrix3d> R_4(parameters_4);
		Eigen::Map<const Eigen::Vector4d> wxyz_4(parameters_4);
		Eigen::Quaterniond q_4(wxyz_4[0], wxyz_4[1], wxyz_4[2], wxyz_4[3]);
		q_4 = q_4.normalized();
		Eigen::Matrix3d R_4 = q_4.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_4(parameters_4 + 4);
		new_node_position += w_j_[4] * (R_4 * (source_ - vector_g_[4]) + vector_g_[4] + t_4);

		//Eigen::Map<const Eigen::Matrix3d> R_5(parameters_5);
		Eigen::Map<const Eigen::Vector4d> wxyz_5(parameters_5);
		Eigen::Quaterniond q_5(wxyz_5[0], wxyz_5[1], wxyz_5[2], wxyz_5[3]);
		q_5 = q_5.normalized();
		Eigen::Matrix3d R_5 = q_5.toRotationMatrix();
		Eigen::Map<const Eigen::Vector3d> t_5(parameters_5 + 4);
		new_node_position += w_j_[5] * (R_5 * (source_ - vector_g_[5]) + vector_g_[5] + t_5);

		res(0, 0) = cost_function_weight_ * (new_node_position - target_)(0);
		res(0, 1) = cost_function_weight_ * (new_node_position - target_)(1);
		res(0, 2) = cost_function_weight_ * (new_node_position - target_)(2);
		return true;
	}

	Eigen::Vector3d source_;
	Eigen::Vector3d target_;
	std::vector<Eigen::Vector3d> vector_g_;
	std::vector<double> w_j_;
	int nodes_connectivity;
	double cost_function_weight_;

};


