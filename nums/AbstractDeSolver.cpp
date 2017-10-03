#include "AbstractDeSolver.h"
#include <cmath>

AbstractDeSolver::AbstractDeSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess)
        : grid_(num_nodes, p_ode->t_min(), p_ode->t_max())
{
    // initialize the ode, grid, and the initial guess
    p_ode_ = p_ode;
    dim_ = p_ode->dim();
    num_nodes_ = num_nodes;
    sol_vec_ = init_guess;

    // file for solution data
    filename_ = "ode_output.txt";
}
//
//void AbstractOdeSolver::SetStepSize(double h)
//{
//    mStepSize = h;
//}
//void AbstractOdeSolver::SetTimeInterval(double t0, double t1)
//{
//    mInitialTime = t0;
//    mFinalTime = t1;
//}
//void AbstractOdeSolver::SetInitialValue(double y0)
//{
//    mInitialValue = y0;
//}
//void AbstractOdeSolver::SetInitialValue(std::vector<double> y0)
//{
//    mInitialValueVector = y0;
//}
//
//double AbstractOdeSolver::time()
//{
//    return t_;
//}
