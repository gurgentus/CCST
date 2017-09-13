//
// Created by Gurgen Hayrapetyan on 9/9/17.
//

#ifndef NUMS_BVPSOLVER_H
#define NUMS_BVPSOLVER_H

#include "Eigen/Dense"
#include "BVPSolver.h"
#include "DifferentialSystem.h"
#include "FiniteDifferenceGrid.h"
#include "Node.h"

class BVPSolver {

private:
    int num_nodes_;  // number of grid nodes
    int dim_;  // dimension of the system
    FiniteDifferenceGrid grid_;  // pointer to a grid

    DifferentialSystem* p_ode_;  // pointer to a system to solve
    Eigen::VectorXd sol_vec_; // pointer to the solution vector
    std::string filename_;  // allow the user to specify the output file or use a default name

    static constexpr double ERR_TOL = 1e-1;
    static const int MAX_ITER = 50;

public:
    BVPSolver(DifferentialSystem *p_ode, int num_nodes, int dim);

    void set_filename(const std::string& name)
    {
        filename_ = name;
    }
    void Solve();
    double Step(int i);
    void WriteSolutionFile();
    Eigen::VectorXd sol_vec() const
    {
        return sol_vec_;
    }
};

#endif //NUMS_BVPSOLVER_H
