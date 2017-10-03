/* author: Gurgen Hayrapetyan.
 * implementation of a Runge-Kutta Solver
 */

#ifndef RUNGEKUTTASOLVER
#define RUNGEKUTTASOLVER

#include "AbstractDeSolver.h"

class RungeKuttaSolver : public AbstractDeSolver {
private:
    double step_size_ = 0.0001;

public:
    RungeKuttaSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess);

    /* returns grid point value */
    double Step(int i) override;

    /* runs the solver */
    int Solve() override;

    void WriteSolutionFile() override;

    Eigen::VectorXd RKIteration(double ti, Eigen::VectorXd yi);
};


#endif //NUMS_RUNGEKUTTASOLVER_H
//    // implementations of virtual methods from inherited class
//    void UpdateState(double dt);
//    void SolveEquation(std::vector<double> yi);
//
//    // single RK iteration
//    void RKIteration(double ti, std::vector<double> &yi);
//    void SetStateDimension(int state_dim);
//
//    void getState(Eigen::VectorXd& st);
//    void setState(const Eigen::VectorXd& st);
//
//    // virtual methods
//    virtual void InitialConditions() = 0;
//    virtual void RightHandSide(double t, const std::vector<double> &  y, std::vector<double> &  f) = 0;

