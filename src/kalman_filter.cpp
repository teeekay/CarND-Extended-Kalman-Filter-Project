#include "kalman_filter.h"
#include <iostream>

using namespace std;

using Eigen::MatrixXd;
using Eigen::VectorXd;

//#define DEBUG
#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif


KalmanFilter::KalmanFilter() {
}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in)
{
	x_ = x_in;
	P_ = P_in;
	F_ = F_in;
	H_ = H_in;
	R_ = R_in;
	Q_ = Q_in;
}

void KalmanFilter::Predict()
{
	x_ = F_ * x_;
	MatrixXd Ft = F_.transpose();
	P_ = F_ * P_ * Ft + Q_;
}

void KalmanFilter::Update(const VectorXd &z)
{
	/**
	   TODO:
	 * update the state by using Kalman Filter equations
	 */

	// KF Measurement Update step
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;

	MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;


	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z)
{
	/**
	   TODO:
	 * update the state by using Extended Kalman Filter equations
	 */
	D(cout << "In UpdateEKF z.size() = " << z.size() << endl; )

	float xcoord = x_(0);
	float ycoord = x_(1);
	float xveloc = x_(2);
	float yveloc = x_(3);

	float rho = sqrt(xcoord * xcoord + ycoord * ycoord);

  /* deal with readings close to zero which can cause divide by zero issues */
	float phi = 0.0;
	if (fabs(xcoord) > 0.0001 )
		phi = atan2(ycoord, xcoord);
	D(cout << "phi = " << phi << endl;)

	float rho_dot = 0.0;
	if(fabs(rho) > 0.0001)
		rho_dot = (xcoord * xveloc + ycoord * yveloc) / rho;
	else
	{
		cout << "adjusting rho to 0.0001 from " << rho << endl;
		rho = 0.0001;
	}

	D(cout << "UpdateEKF: set up z_pred" << endl;)
	VectorXd z_pred = VectorXd(3);
	z_pred << rho, phi, rho_dot;

	D(cout << "UpdateEKF: z_pred = " << z_pred << endl;)

	VectorXd y = z - z_pred;

  /* make sure that y(1) is between -PI and PI */
  while(y(1) > M_PI || y(1) < -M_PI)
	{
		cout << "UpdateEKF y(1) out of range = " << y(1) << endl;
		while(y(1) > M_PI)
			y(1) -= 2 * M_PI;
		while(y(1) < -M_PI)
			y(1) += 2 * M_PI;
	}

	D(cout << "UpdateEKF: y =" << y << endl;)

  MatrixXd Ht = H_.transpose();
	MatrixXd PHt = P_ * Ht;
	MatrixXd S = H_ * PHt + R_;
	MatrixXd Si = S.inverse();
	MatrixXd K = PHt * Si;

	x_ = x_ + (K * y);
	P_ = (I_ - K * H_) * P_;
}
