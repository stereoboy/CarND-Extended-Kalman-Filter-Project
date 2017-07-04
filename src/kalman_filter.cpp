#include "kalman_filter.h"
#include <cmath>
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

#define MIN(X,Y) ((X) < (Y) ? (X) : (Y))  

KalmanFilter::KalmanFilter() {}

KalmanFilter::~KalmanFilter() {}

void KalmanFilter::Init(VectorXd &x_in, MatrixXd &P_in, MatrixXd &F_in,
                        MatrixXd &H_in, MatrixXd &R_in, MatrixXd &Q_in) {
  x_ = x_in;
  P_ = P_in;
  F_ = F_in;
  H_ = H_in;
  R_ = R_in;
  Q_ = Q_in;
}

void KalmanFilter::Predict() {
  /**
  TODO:
    * predict the state
  */
  x_ = F_ * x_;
  MatrixXd Ft = F_.transpose();
  P_ = F_ * P_ * Ft + Q_; 
}

void KalmanFilter::Update(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Kalman Filter equations
  */
	VectorXd z_pred = H_ * x_;
	VectorXd y = z - z_pred;
	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}

void KalmanFilter::UpdateEKF(const VectorXd &z) {
  /**
  TODO:
    * update the state by using Extended Kalman Filter equations
  */
  float px = x_[0];
  float py = x_[1];
  float vx = x_[2];
  float vy = x_[3];
  
	VectorXd z_pred = VectorXd(3);
  z_pred << sqrt(px*px + py*py), atan2(py, px), (px*vx + py*vy)/sqrt(px*px + py*py);

	VectorXd y = z - z_pred;
  //normalize angle difference
  cout<< "y:\n" << y << endl;
  if (y[1] < -M_PI)
  {
    //cout<<"CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC" << endl;
    y[1] += 2*M_PI;
  }
  else if (y[1] > M_PI)
  {
    //cout<<"DDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD" << endl;
    y[1] -= 2*M_PI;
  }

#if 0
  cout<< "x:\n" << x_ << endl;
  cout<< "z vs z_pred\n";
  cout<< z << "\n" << "pred:\n" << z_pred << endl;
  cout<< "y:\n" << y << endl;

  cout<< "px,py=" << px <<", "<< py << endl;
  cout<< "atan2(py, px)=" << atan2(py, px) << endl;
  cout<<"======================"<<endl;
#endif

	MatrixXd Ht = H_.transpose();
	MatrixXd S = H_ * P_ * Ht + R_;
	MatrixXd Si = S.inverse();
	MatrixXd PHt = P_ * Ht;
	MatrixXd K = PHt * Si;

	//new estimate
	x_ = x_ + (K * y);
	long x_size = x_.size();
	MatrixXd I = MatrixXd::Identity(x_size, x_size);
	P_ = (I - K * H_) * P_;
}
