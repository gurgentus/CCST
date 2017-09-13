#include <iostream>
#include "DifferentialSystem.h"
#include "BVPSolver.h"
#include "OrbitTransfer.h"

Eigen::VectorXd rhs(double t, const Eigen::VectorXd& y)
{
//    Eigen::VectorXd v = Eigen::VectorXd(2);
//    v << y(1), 34*sin(t)+4*y(0)-3*y(1);
    Eigen::VectorXd rs = Eigen::VectorXd(7);
    double r = y(0);
    double u = y(1);
    double v = y(2);

    double pr = y(4);
    double pu = y(5);
    double pv = y(6);
    //double pth = y(7);
    double rsq = r*r;
    double m0 = 400;
    double tf = 4*24*3600;
    double mdottau = 0.5*m0/tf;
    double r0 = 6678;
    double v0 = 7.726;
    double eta = v0*tf/r0;
    double T = 20*tf/v0;
    double post = (T/(m0-mdottau*t*tf))/(sqrt(pu*pu+pv*pv));

    rs << u*eta,
            (v*v/r-1/rsq)*eta+post*pu,
            -u*v*eta/r + post*pv,
            v*eta/r,
            pu*(v*v/rsq-2/(rsq*r))*eta-pv*u*v*eta/rsq,//+pth*v*eta/rsq,
            -pr*eta+pv*v/r,
            -pu*2*v*eta/r+pv*u*eta/r;//-pth*eta/r;
    //0;
    return rs;
}

Eigen::MatrixXd rhs_grad(double t, const Eigen::VectorXd& y)
{
//    Eigen::MatrixXd v = Eigen::MatrixXd(2,2);
//    v << 0,1,
//         4,-3;

    double r = y(0);
    double u = y(1);
    double v = y(2);

    double pu = y(5);
    double pv = y(6);
    //double pth = y(7);
    double rsq = r*r;
    double m0 = 400;
    double tf = 4*24*3600;
    double mdottau = 0.5*m0/tf;
    double r0 = 6678;
    double v0 = 7.726;
    double eta = v0*tf/r0;
    double T = 20*tf/v0;
    double post = (T/(m0-mdottau*t*tf))/(sqrt(pu*pu+pv*pv));

    double pupvsq = pu*pu+pv*pv;
    double prdr = pu*(-2*v*v/(rsq*r)+6/(rsq*rsq))*eta +2*pv*u*v*eta/(rsq*r);//-pth*2*v*eta/(r*rsq);
    double prdu = -pv*v*eta/rsq;
    double prdv =  2*pu*v*eta/rsq-pv*u*eta/rsq;//+pth*eta/rsq;

    Eigen::MatrixXd rs = Eigen::MatrixXd(7,7);
    rs << 0,eta,0,0,0,0,0,//0,
            eta*(-v*v/rsq+2/(rsq*r)), 0, 2*v*eta/r, 0, 0, post-post*pu*pu/pupvsq, -post*pu*pv/pupvsq,//0,
            u*v*eta/rsq, -v*eta/r, -u*eta/r, 0, 0, -post*pu*pv/pupvsq, post-post*pv*pv/pupvsq,//0,
            -v*eta/rsq, 0, eta/r, 0,0,0,0,//0,
            prdr, prdu, prdv, 0, 0, (v*v/rsq-2/(rsq*r))*eta, -u*v*eta/rsq, //v*eta/rsq,
            -pv*v/rsq, 0, pv/r, 0, -eta, 0, v/r, //0,
            //pu*2*v*eta/rsq-pv*u*eta/rsq+pth*eta/rsq, pv*eta/r, -pu*2*eta/r, 0, 0, -2*v*eta/r, u*eta/r, -eta/r,
            pu*2*v*eta/rsq-pv*u*eta/rsq, pv*eta/r, -pu*2*eta/r, 0, 0, -2*v*eta/r, u*eta/r;
    //0,0,0,0,0,0,0,0;

    return rs;
}
Eigen::VectorXd bc(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
    double r = y1(0);
    double u = y1(1);
    double v = y1(2);
    double theta = y1(3);

    double r2 = y2(0);
    double u2 = y2(1);
    double v2 = y2(2);
    double pr2 = y2(4);

    double pv2 = y2(6);
    //double pth2 = y2(7);
//    Eigen::VectorXd v = Eigen::VectorXd(2);
//    v << y1(1)+5, y2(0)-4;
    Eigen::VectorXd rs = Eigen::VectorXd(7);
    rs << r-1, u, v-1, theta, u2, v2-sqrt(1/r2), 1-pr2+0.5*pv2*sqrt(1/(r2*r2*r2));//, pth2;
    return rs;
}

Eigen::MatrixXd bc_grad1(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
//    Eigen::MatrixXd v = Eigen::MatrixXd(2,2);
//    v << 0,1,
//         0,0;
    Eigen::MatrixXd rs = Eigen::MatrixXd(7,7);
    rs << 1,0,0,0,0,0,0,//0,
            0,1,0,0,0,0,0,//0,
            0,0,1,0,0,0,0,//0,
            0,0,0,1,0,0,0,//0,
            0,0,0,0,0,0,0,//0,
            0,0,0,0,0,0,0,//0,
            0,0,0,0,0,0,0;//0,
    //0,0,0,0,0,0,0,0;

    return rs;
}

Eigen::MatrixXd bc_grad2(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
//    Eigen::MatrixXd v = Eigen::MatrixXd(2,2);
//    v << 0,0,
//         1,0;

    double r2 = y2(0);
    double pv2 = y2(6);

    Eigen::MatrixXd rs = Eigen::MatrixXd(7,7);
//    rs << 0,0,0,0,0,0,0,0,
//          0,0,0,0,0,0,0,0,
//          0,0,0,0,0,0,0,0,
//          0,0,0,0,0,0,0,0,
//          0,1,0,0,0,0,0,0,
//          0.5*pow(r2,-1.5),0,1,0,0,0,0,0,
//          -(3.0/4.0)*pv2*pow(r2,-2.5),0,0,0,-1,0,0.5*sqrt(1/(r2*r2*r2)),0,
//          0,0,0,0,0,0,0,1;
    rs << 0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,
            0,0,0,0,0,0,0,
            0,1,0,0,0,0,0,
            0.5*pow(r2,-1.5),0,1,0,0,0,0,
            -(3.0/4.0)*pv2*pow(r2,-2.5),0,0,0,-1,0,0.5*sqrt(1/(r2*r2*r2));
    return rs;
}

int main() {
//    std::cout << "Hello, World!" << std::endl;
//
//    std::cout << "Setting up DS" << std::endl;
//    DifferentialSystem ds(rhs, rhs_grad, bc, bc_grad1, bc_grad2, 0.0, 1.0);
//    std::cout << "Setting up BVP" << std::endl;
//    BVPSolver bvp(&ds, 101, 7);
//    std::cout << "Solving..." << std::endl;
//    bvp.Solve();
//    std::cout << "Solved." << std::endl;
//    for (int i=0; i < bvp.sol_vec().size()/2; i++)
//    {
//        std::cout << bvp.Step(i) << bvp.sol_vec()(2*i) << std::endl;
//    }

    OrbitTransfer ot;
    ot.run();

    return 0;
}

