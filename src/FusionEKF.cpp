#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>
//#define DEBUG

#ifdef DEBUG
#define D(x) x
#else
#define D(x)
#endif

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF()
{
  is_initialized_ = false;
  previous_timestamp_ = 0;

  // initializing matrices
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
    0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
    0, 0.0009, 0,
    0, 0, 0.09;

  /**
     TODO:
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
  H_laser_ << 1, 0, 0, 0,   /*from section 10 in lesson 5 */
    0, 1, 0, 0;

  ekf_.F_ = MatrixXd(4, 4); /* 4x4 matrix (state transition) */

  ekf_.P_ = MatrixXd(4, 4); /* 4x4 matrix  covariamce */

  ekf_.I_ = MatrixXd(4, 4); /* 4x4 Identity matrix */

  /* set the acceleration noise components */
  noise_ax = 9.0; /* set to 5 in lesson 5 section 13 - use 5 as directed below*/
  noise_ay = 9.0; /* set to 5 in lesson 5 section 13 - use 5 as directed below*/
  /**
     TODO:
   * Finish initializing the FusionEKF.
   * Set the process and measurement noises
   */
}

/**
 * Destructor.
 */
FusionEKF::~FusionEKF()
{
}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack)
{
  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_)
    {
      /**
         TODO:
       * Initialize the state ekf_.x_ with the first measurement.
       * Create the covariance matrix.
       * Remember: you'll need to convert radar from polar to cartesian coordinates.
       */
      // first measurement
      cout << "FusionEKF Initialization " << endl;
      ekf_.x_ = VectorXd(4);
      ekf_.x_ << 1, 1, 1, 1; /* play with these they will affect RMSE */
      D(cout << "FusionEKF set x_" << endl;)

      if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
        {
          D(cout << "FusionEKF Radar Measurement" << endl;)
          /**
             Convert radar from polar to cartesian coordinates and initialize state.
           */
          // just set ekf_.x_(0) to ro*cos(theta)
          ekf_.x_(0) = measurement_pack.raw_measurements_(0) * cos(measurement_pack.raw_measurements_(1));
          // just set ekf_.x_(1) to ro*sin(theta)
          ekf_.x_(1) = measurement_pack.raw_measurements_(0) * sin(measurement_pack.raw_measurements_(1));
        }
      else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER)
        {
          D(cout << "FusionEKF Laser Measurement" << endl;)
          /**
             Initialize state.
           */
          // set ekf_.x_(0)to x
          ekf_.x_(0) = measurement_pack.raw_measurements_(0);
          // set ekf_.x_(1)to y
          ekf_.x_(1) = measurement_pack.raw_measurements_(1);
        }
      //ekf_.F_ << set to 1 diagonal matrix
      D(cout <<"initialize F_"<<endl;)
      ekf_.F_ << 1, 0, 1, 0,
                 0, 1, 0, 1,
                 0, 0, 1, 0,
                0, 0, 0, 1;

      D(cout <<"initialize P_"<<endl;)
      ekf_.P_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1000, 0,
                 0, 0, 0, 1000;
      D(cout <<"initialize I_"<<endl;)
      ekf_.I_ << 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1;

      previous_timestamp_ = measurement_pack.timestamp_;

      // done initializing, no need to predict or update
      is_initialized_ = true;
      cout << "FusionEKF Initialized - timestamp: " << previous_timestamp_ << endl;
      return;
    }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     TODO:
   * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
   * Update the process noise covariance matrix.
   * Use noise_ax = 9 and noise_ay = 9 for your Q matrix. -> (TK set above)
   */
  D(cout << "FusionEKF - calc dt" << endl;)
  float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;
  previous_timestamp_ = measurement_pack.timestamp_;

  float dt_2 = dt * dt;
  float dt_3 = dt_2 * dt;
  float dt_4 = dt_3 * dt;

  D(cout << "FusionEKF - modify F_" << endl;)
  //Modify the F matrix so that the time is integrated
  ekf_.F_(0, 2) = dt;
  ekf_.F_(1, 3) = dt;

  D(cout << "FusionEKF - set up Q_" << endl;)
  //set the process covariance matrix 0 section 9 of lesson 5
  ekf_.Q_ = MatrixXd(4, 4);
  ekf_.Q_ << dt_4 / 4 * noise_ax, 0, dt_3 / 2 * noise_ax, 0,
    0, dt_4 / 4 * noise_ay, 0, dt_3 / 2 * noise_ay,
    dt_3 / 2 * noise_ax, 0, dt_2 * noise_ax, 0,
    0, dt_3 / 2 * noise_ay, 0, dt_2 * noise_ay;

  D(cout << "FusionEKF - call Predict()" << endl;)
  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     TODO:
   * Use the sensor type to perform the update step.
   * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR)
    {
      cout << "FusionEKF update with Radar Measurement" << endl;
      // Radar updates
      //set ekf_.H_ by setting to Hj which is the calculated jacobian (look in tools.h)
      ekf_.H_ = tools.CalculateJacobian(ekf_.x_);
      ekf_.R_ = R_radar_;
      D(cout << "FusionEKF call UpdateEKF with Radar Measurement" << endl;)
      ekf_.UpdateEKF(measurement_pack.raw_measurements_);
    }
  else
    {
      // Laser updates
      cout << "FusionEKF update with Laser Measurement" << endl;
      ekf_.H_ = H_laser_;
      ekf_.R_ = R_laser_;
      D(cout << "FusionEKF call Update with Laser Measurement" << endl;)
      ekf_.Update(measurement_pack.raw_measurements_);
    }

  // print the output
  cout << "x_ = " << ekf_.x_ << endl;
  cout << "P_ = " << ekf_.P_ << endl;
}
