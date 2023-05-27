/* author: Gurgen Hayrapetyan.
 * implementation of a collocation with quasilinearization
 * ref: U.Ascher, R. Mattheij, R. Russelll, Numerical Solution of Boundary Value Problems
 * for Ordinary Differential Equations
 * The scheme defaults to a three stage Lobatto method, but can be reset using the SetScheme
 * method.
 * See OrbitTransfer.cpp for usage example.
 */


#ifndef NUMS_COLLOCATIONSOLVER_H
#define NUMS_COLLOCATIONSOLVER_H

#include "AbstractDeSolver.h"

class CollocationSolver : public AbstractDeSolver {
private:
    /* Collocation Parameters */
    int k = 3;
    Eigen::MatrixXd a = Eigen::MatrixXd(3,3);
    Eigen::VectorXd rho = Eigen::VectorXd(3);
    Eigen::VectorXd b = Eigen::VectorXd(3);

public:
    CollocationSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess);
    virtual ~CollocationSolver() {};
    
    /* sets the collocation scheme */
    void SetScheme(int k, Eigen::VectorXd& rho, Eigen::MatrixXd& a, Eigen::VectorXd& b);

    /* returns grid point value */
    double Step(int i) override;

    /* runs the solver */
    int Solve() override;

    void WriteSolutionFile() override;
};


#endif //NUMS_COLLOCATIONSOLVER_H
