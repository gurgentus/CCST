#include <iostream>
#include "DifferentialSystem.h"
#include "BVPSolver.h"
#include "Car.h"


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

//    OrbitTransfer ot;
//    ot.run();

    Eigen::VectorXd init_state = Eigen::VectorXd::Zero(5);
    init_state << -100, 0, 5, 0, 1;

    Car car(5, 500);

    std::vector<Eigen::MatrixXd> waypoints(2);
//    waypoints[0] = Eigen::MatrixXd(3,3);
//    waypoints[0] << -100, -50, 0,
//                    5, 5, 5,
//                    0, 0, 0;
//    waypoints[1] = Eigen::MatrixXd(3,3);
//    waypoints[1] << 0, 4, 0,
//                    0, 0, 0,
//                    0, 0, 0;
//
//    Eigen::VectorXd time_at_waypoint = Eigen::VectorXd(3);
//    time_at_waypoint << 0, 10, 20;

    waypoints[0] = Eigen::MatrixXd(3,2);
    waypoints[0] << 909.654, 970.647,
            13.20148, 1.94931,
            0, 0;
    waypoints[1] = Eigen::MatrixXd(3,2);
    waypoints[1] << 1128.67, 1133.51,
            0, 0.447428,
            0, 0;


    Eigen::VectorXd time_at_waypoint = Eigen::VectorXd(2);
    time_at_waypoint << 0, 10;

    int status = car.SimulateMinJerkLaneChange(init_state, waypoints, time_at_waypoint);
    if (status == 2) {
        std::cout << "Solver failed" << std::endl;
    }
    return 0;
}

