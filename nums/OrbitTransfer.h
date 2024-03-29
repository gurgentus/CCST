/* author: Gurgen Hayrapetyan
 * maximum radius orbit transfer optimal control problem
 * ref: James M Longuski, José J. Guzmán, John E. Prussing (auth.)
 *  -Optimal Control with Aerospace Applications-Springer-Verlag New York (2014)
 */


#ifndef NUMS_ORBITTRANSFER_H
#define NUMS_ORBITTRANSFER_H

#include "Eigen/Dense"
#include "DifferentialSystem.h"
#include <boost/python/numpy.hpp>

class OrbitTransfer : public DifferentialSystem {
private:
    double mu = 3.986004418e5; // gravitational constant
    double m0 = 1500; // mass (kg)
    double r0 = 6378+300; // initial radius
    double v0 = sqrt(mu/r0); // initial tangential velocity of circular orbit
    double T = 20; // thrust (N)
    double tf = 0.1*24*3600; // final time (s)
    double Tbar = T*tf/(v0*1000); // nondimensionalized thrust
    double Isp = 6000; // specific impulse
    double mdot = T/Isp/9.80665; // mass flow rate
    double eta = v0*tf/r0; // scaling

    int N = 500;
    int k = 3;
    Eigen::MatrixXd a = Eigen::MatrixXd(3,3);
    Eigen::VectorXd rho = Eigen::VectorXd(3);
    Eigen::VectorXd b = Eigen::VectorXd(3);
    Eigen::VectorXd sol_vec_;
public:
    OrbitTransfer();
    OrbitTransfer(double mu, double m0, double Isp, double T, double r0, double tf, int N);
    Eigen::VectorXd RhsFunc(double t, const Eigen::VectorXd& y) override;
    Eigen::MatrixXd RhsGradYFunc(double t, const Eigen::VectorXd& y) override;
    Eigen::MatrixXd BcsGrad1Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;
    Eigen::MatrixXd BcsGrad2Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;
    Eigen::VectorXd BcsFunc(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2) override;
    int run(double mu, double m0, double Isp, double T, double r0, double days, double timestep_hrs, int N);
    int SetMatrix(int k, boost::python::list& rho, boost::python::list& a, boost::python::list& b);
    boost::python::list getT();
    boost::python::list getAngle();
    boost::python::list getX();
    boost::python::list getY();
    char const* greet();
};


#endif //NUMS_ORBITTRANSFER_H
