#include "tools.h"
#include <iostream>
//#define DEBUG

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools()
{
}

Tools::~Tools()
{
}


VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth)
{
	VectorXd rmse(4);

	rmse << 0, 0, 0, 0;
	//  check the validity of the following inputs:
	//  the estimation vector size should not be zero
	//  * the estimation vector size should equal ground truth vector size
	if ((estimations.size() != ground_truth.size()) ||
	    (estimations.size() == 0))
	{
		cout << "Error - CalculateRMSE: estimations.size() =" <<
		estimations.size()
		     << "ground_truth.size() = " << ground_truth.size() << endl;
		return rmse;
	}

	//accumulate squared residuals
	for (int i = 0; i < estimations.size(); ++i)
	{
		VectorXd residual = estimations[i] - ground_truth[i];
		residual = residual.array() * residual.array();
		rmse += residual;
	}
	//calculate the mean
	rmse = rmse / estimations.size();
	//calculate the squared root
	rmse = sqrt(rmse.array());
	//return the result
	return rmse;
}

MatrixXd Tools::CalculateJacobian(const VectorXd& x_state)
{
	/**
	   TODO:
	 * Calculate a Jacobian here.
	 */
	MatrixXd Hj(3, 4);

	Hj << 0, 0, 0, 0,
	        0, 0, 0, 0,
	        0, 0, 0, 0;
	//recover state parameters
	float px = x_state(0);
	float py = x_state(1);
	float vx = x_state(2);
	float vy = x_state(3);

	float pxy_sqrd = px * px + py * py;

	//check division by zero
	if (pxy_sqrd <= 0.0001)
	{
		cout << "Error CalculateJacobian: divide by zero" << endl;
		return Hj;
	}

	float pxy_sqrt = sqrt(pxy_sqrd);
	float pxy_sqrtd = pxy_sqrd * pxy_sqrt;
	Hj << px / pxy_sqrt, py / pxy_sqrt, 0, 0,
	        -py / pxy_sqrd, px / pxy_sqrd, 0, 0,
	        py * (vx * py - vy * px) / pxy_sqrtd,
	        px * (vy * px - vx * py) / pxy_sqrtd, px / pxy_sqrt, py / pxy_sqrt;

	D(cout << "tools: Hj = " << Hj << endl;)

	return Hj;
}
