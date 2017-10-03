/* author: Gurgen Hayrapetyan
 * car dynamics model and trajectory planners
 */

#ifndef ORBITTRANSFER_CAR_H
#define ORBITTRANSFER_CAR_H

#include "Eigen/Dense"
#include "DifferentialSystem.h"
#include <vector>
//#include <boost/numpy.hpp>

class Car : public DifferentialSystem {
private:
    int N = 500; // number of gridpoints for solvers and simulations
    double l = 5; // distance between front and back axes
    Eigen::VectorXd sol_vec_; // saved trajectory

    // current state for simulations
    double x = 0;
    double y = 0;
    double v = 0;
    double xi = 0;
    double m = 0;
    // current control values for simulations
    double dw_ = 0;
    double T_ = 0;

    enum ControlScheme {MinTime, MinJerk} control_scheme;

public:
    Car(double l, int N);
    /* Dynamics */
    Eigen::VectorXd RhsFunc(double t, const Eigen::VectorXd& y) override;
    Eigen::MatrixXd RhsGradYFunc(double t, const Eigen::VectorXd& y) override;
    Eigen::MatrixXd BcsGrad1Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;
    Eigen::MatrixXd BcsGrad2Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;
    Eigen::VectorXd BcsFunc(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;

    /* Control schemes */
    void SetControlScheme(ControlScheme scheme);
    void UpdateControls(const Eigen::VectorXd& y_vec, ControlScheme control_scheme);
    double MinTimeStoppageThrottleControl(const Eigen::VectorXd& y_vec);
    int GenerateMinDerivTraj(const Eigen::MatrixXd& waypoints, Eigen::VectorXd& spline_coeffs);
    int GenerateMinDerivTraj(int num_coeffs, std::vector<Eigen::MatrixXd>& waypoints, Eigen::VectorXd& time_at_waypoint, Eigen::MatrixXd& spline_coeffs);
    void DesiredState(int num_coeffs, double t, Eigen::MatrixXd& spline_coeffs, const Eigen::VectorXd& time_at_waypoint, Eigen::VectorXd& state);
    int SimulateMinJerkLaneChange(const Eigen::VectorXd& init_y, std::vector<Eigen::MatrixXd>& waypoints, Eigen::VectorXd& time_at_waypoint);
    int SimulateMinTimeBreaking(const Eigen::VectorXd& init_y);
    int Factorial(int n);
    boost::python::list Generate2DTrajectory(boost::python::list& start_pos, boost::python::list& start_vel,
                                                  boost::python::list& start_acc, boost::python::list& end_pos,
                                                  boost::python::list& end_vel, boost::python::list& end_acc, double time_ahead);
//    int run(double mu, double m0, double Isp, double T, double r0, double days, double timestep_hrs, int N);
//    int SetMatrix(int k, boost::python::list& rho, boost::python::list& a, boost::python::list& b);
//    boost::python::list getT();
//    boost::python::list getAngle();
//    boost::python::list getX();
//    boost::python::list getY();
//    char const* greet();
};

#endif //ORBITTRANSFER_CAR_H
