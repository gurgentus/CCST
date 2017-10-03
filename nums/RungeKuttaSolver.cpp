/* author: Gurgen Hayrapetyan.
 * implementation of a Runge-Kutta Solver
 */

#include "RungeKuttaSolver.h"
#include <iostream>
#include <fstream>
#include <cassert>
#include <cmath>

RungeKuttaSolver::RungeKuttaSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess)
        : AbstractDeSolver(p_ode, num_nodes, init_guess)
{
    // initialize the ode, grid, and the initial guess
    p_ode_ = p_ode;
    dim_ = p_ode->dim();
    num_nodes_ = num_nodes;
    int k = 1;
    sol_vec_ = Eigen::VectorXd((num_nodes_+1)*k*dim_);
    sol_vec_.segment(0,dim_) = init_guess;

    // file for solution data
    filename_ = "ode_output.txt";

    f_vec = Eigen::VectorXd::Zero((num_nodes_+1)*k*dim_);
    for (int i=0; i<num_nodes_; i++)
    {
        f_vec.block(i*k*dim_, 0, dim_, 1) = p_ode_->RhsFunc(grid_.nodes_[i].coordinate, sol_vec_.segment(i*dim_, dim_));
    }

    step_size_ = (p_ode->t_max() - p_ode->t_min())/num_nodes_;
}

double RungeKuttaSolver::Step(int i)
{
    return grid_.nodes_[i].coordinate;
}


Eigen::VectorXd RungeKuttaSolver::RKIteration(double ti, Eigen::VectorXd yi)
{
    Eigen::VectorXd k1 = Eigen::VectorXd(dim_);
    Eigen::VectorXd k2 = Eigen::VectorXd(dim_);
    Eigen::VectorXd k3 = Eigen::VectorXd(dim_);
    Eigen::VectorXd k4 = Eigen::VectorXd(dim_);
    Eigen::VectorXd f1 = Eigen::VectorXd(dim_);
    Eigen::VectorXd f2 = Eigen::VectorXd(dim_);
    Eigen::VectorXd f3 = Eigen::VectorXd(dim_);
    Eigen::VectorXd f4 = Eigen::VectorXd(dim_);
    Eigen::VectorXd int2 = Eigen::VectorXd(dim_);
    Eigen::VectorXd int3 = Eigen::VectorXd(dim_);
    Eigen::VectorXd int4 = Eigen::VectorXd(dim_);

    //RightHandSide(ti, yi, f1);
    f1 = p_ode_->RhsFunc(ti, yi);
    for (int j=0; j < dim_; j++)
    {
        k1(j) = step_size_*f1(j);
        int2(j) = yi(j) + 0.5 * k1(j);
    }
    //RightHandSide(ti+0.5*step_size_, int2, f2);
    f2 = p_ode_->RhsFunc(ti+0.5*step_size_, int2);
    for (int j=0; j < dim_; j++)
    {
        k2(j) = step_size_*f2(j);
        int3(j) = yi(j) + 0.5 * k2(j);
    }
    //RightHandSide(ti+0.5*step_size_, int3, f3);
    f3 = p_ode_->RhsFunc(ti+0.5*step_size_, int3);
    for (int j=0; j < dim_; j++)
    {
        k3(j) = step_size_*f3[j];
        int4(j) = yi(j)+k3(j);
    }

    //RightHandSide(ti+step_size_, int4, f4);
    f4 = p_ode_->RhsFunc(ti+step_size_, int4);
    for (int j=0; j < dim_; j++)
    {
        k4(j) = step_size_*f4(j);
        yi(j) = yi(j) + (1.0/6.0)*(k1(j) + 2*k2(j) + 2*k3(j) + k4(j));
    }
    return yi;

}

int RungeKuttaSolver::Solve()
{
//    double numPoints = (mFinalTime - mInitialTime)/mStepSize;
//    double ti = mInitialTime;
    double ti = p_ode_->t_min();
    for (int i=0; i < num_nodes_; i++) {
        sol_vec_.segment((i + 1) * dim_, dim_) = RKIteration(ti, sol_vec_.segment(i*dim_, dim_));
        ti = ti + step_size_;
    }
    return 0;
}

void RungeKuttaSolver::WriteSolutionFile()
{
    std::ofstream output_file(filename_.c_str());
    assert(output_file.is_open());
    for (int i=0; i<num_nodes_; i++)
    {
        double t = grid_.nodes_[i].coordinate;
        output_file << t << " " << sol_vec_(i*dim_) << "\n";
    }
    output_file.flush();
    output_file.close();
    std::cout << "Solution written to " << filename_ << std::endl;

}
//
//void RungeKuttaSolver::UpdateState(double dt)
//{
//    double total_d = 0;
//    while (total_d < dt)
//    {
//        RKIteration(t_, state);
//        total_d += mStepSize;
//    }
//    t_ = t_ + total_d;
//}
//
//void RungeKuttaSolver::getState(Eigen::VectorXd& st)
//{
//    st = Eigen::VectorXd(state_dim_);
//    for (int i=0; i<state_dim_; i++ )
//    {
//         st(i) = state[i];
//    }
//}
//
//void RungeKuttaSolver::setState(const Eigen::VectorXd& st)
//{
//    for (int i=0; i<state_dim_; i++ )
//    {
//         state[i] = st(i);
//    }
//}
//
//
//
//
