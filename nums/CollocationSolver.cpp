/* author: Gurgen Hayrapetyan.
 * implementation of a collocation with quasilinearization
 * ref: U.Ascher, R. Mattheij, R. Russelll, Numerical Solution of Boundary Value Problems
 * for Ordinary Differential Equations
 */

#include "CollocationSolver.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include "Eigen/SparseLU"

CollocationSolver::CollocationSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess)
        : grid_(num_nodes, p_ode->t_min(), p_ode->t_max())
{
    // initialize the ode, grid, and the initial guess
    p_ode_ = p_ode;
    dim_ = p_ode->dim();
    num_nodes_ = num_nodes;
    sol_vec_ = init_guess;

    // file for solution data
    filename_ = "ode_output.txt";

    // default collocation parameters - three stage Lobatto scheme
    k = 3;
    rho << 0.0, 0.5, 1.0;
    b << 1.0/6.0, 2.0/3.0, 1.0/6.0;
    a << 0, 0, 0,
            5.0/24.0, 1.0/3.0, -1.0/24.0,
            1.0/6.0, 2.0/3.0, 1.0/6.0;

    f_vec = Eigen::VectorXd::Zero((num_nodes_+1)*k*dim_);
    for (int i=0; i<num_nodes_; i++)
    {
        //for (int j=0; j<k; j++) {
        //    f_vec.block(i * k * dim_ + j*dim_, 0, dim_, 1) = p_ode_->RhsFunc(grid_.nodes_[i].coordinate,
        //                                                            sol_vec_.segment(i * dim_, dim_));
        //}
        f_vec.block(i*k*dim_, 0, dim_, 1) = p_ode_->RhsFunc(grid_.nodes_[i].coordinate, sol_vec_.segment(i*dim_, dim_));
    }
}

void CollocationSolver::SetScheme(int k, Eigen::VectorXd& rho, Eigen::MatrixXd& a, Eigen::VectorXd& b) {
    this->k = k;
    this->rho = rho;
    this->a = a;
    this->b = b;
}

double CollocationSolver::Step(int i)
{
    return grid_.nodes_[i].coordinate;
}

int CollocationSolver::Solve()
{
    for (int iter=0; iter<MAX_ITER; iter++)
    {
        std::cout << "Iteration: " << iter << std::endl;

        // Setup matrices involved in the method
        Eigen::VectorXd u = Eigen::VectorXd(sol_vec_.segment(0,dim_));
        Eigen::VectorXd v = Eigen::VectorXd(sol_vec_.tail(dim_));

        Eigen::MatrixXd B1 = p_ode_->BcsGrad1Func(u,v);
        Eigen::MatrixXd B2 = p_ode_->BcsGrad2Func(u,v);
        Eigen::VectorXd beta = -p_ode_->BcsFunc(u,v);

        Eigen::MatrixXd A = Eigen::MatrixXd(dim_, dim_);
        Eigen::MatrixXd W = Eigen::MatrixXd(k*dim_, k*dim_);
        Eigen::MatrixXd V = Eigen::MatrixXd(k*dim_, dim_);
        Eigen::VectorXd q = Eigen::VectorXd(k*dim_);
        Eigen::VectorXd r = Eigen::VectorXd((num_nodes_+1)*dim_);
        Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim_,dim_);
        Eigen::MatrixXd I2 = Eigen::MatrixXd::Identity(k*dim_,k*dim_);
        Eigen::MatrixXd D = Eigen::MatrixXd(dim_, dim_*k);

        Eigen::MatrixXd U = Eigen::MatrixXd(k*num_nodes_*dim_, dim_);
        Eigen::VectorXd p = Eigen::VectorXd(k*num_nodes_*dim_);

        for (int i=0; i<k; i++)
        {
            D.block(0, i*dim_, dim_, dim_) = b(i)*I;
        }

        // This will be a large sparse matrix so use Eigen::SparseMatrix
        Eigen::SparseMatrix<double, Eigen::ColMajor> Gamma((num_nodes_+1)*dim_, (num_nodes_+1)*dim_);

        for (int i=0; i<num_nodes_; i++)
        {
            double t_i = grid_.nodes_[i].coordinate;
            double t_ip1 = grid_.nodes_[i+1].coordinate;
            double h = t_ip1-t_i;
            double hi = h/(k-1);

            Eigen::VectorXd t = Eigen::VectorXd(k);
            Eigen::VectorXd yi = Eigen::VectorXd(sol_vec_.segment(i*dim_,dim_));

            std::vector<Eigen::VectorXd> y;

            for (int j=0; j<k; j++)
            {
		t(j) = t_i + hi*rho(j);
                Eigen::VectorXd yi_temp = yi;

                for (int l=0; l<k; l++)
                {
                    yi_temp = yi_temp + hi*a(j,l)*f_vec.segment((i*k+l)*dim_, dim_);
                }

                A = p_ode_->RhsGradYFunc(t(j),yi_temp);

                for (int l=0; l<k; l++)
                {
                    W.block(j*dim_,l*dim_, dim_, dim_) = a(j,l)*A;
                }

                V.block(j*dim_, 0, dim_, dim_) = A;
                q.block(j*dim_, 0, dim_, 1) = p_ode_->RhsFunc(t(j), yi_temp) - f_vec.segment((i*k+j)*dim_, dim_);
            }


            W = I2 - hi*W;

            Eigen::MatrixXd Winv;
            Eigen::MatrixXd S;
            Winv = W.inverse();
            U.block(i*k*dim_, 0, k*dim_, dim_) = Winv*V;
            p.segment(i*k*dim_, k*dim_) = Winv*q;
            S = -I - hi*D*Winv*V;

            // populate the sparse matrix
            for (int i2=0; i2<dim_; i2++)
            {
                for (int j2=0; j2<dim_; j2++)
                {
                    Gamma.coeffRef(i*dim_+i2, i*dim_+j2) = S(i2,j2);
                    Gamma.coeffRef(i*dim_+i2, (i+1)*dim_+j2) = I(i2,j2);
                }
            }
            r.block(i*dim_, 0, dim_, 1) = hi*D*Winv*q;
        }

        // populate the sparse matrix
        for (int i2=0; i2<dim_; i2++)
        {
            for (int j2=0; j2<dim_; j2++)
            {

                Gamma.coeffRef(num_nodes_*dim_+i2, j2) = B1(i2,j2);
                Gamma.coeffRef(num_nodes_*dim_+i2, num_nodes_*dim_+j2) = B2(i2,j2);
            }
        }
        r.block(num_nodes_*dim_, 0, dim_, 1) = beta;

        // prepare for solving the system with sparse matrix
        Gamma.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        // Compute the ordering permutation vector from the structural pattern of Gamma
        solver.analyzePattern(Gamma);
        // Compute the numerical factorization
        solver.factorize(Gamma);
        //solver.compute(Gamma);

        std::cout << "Solver info: " << solver.info() << std::endl;
	if (solver.info()  == 1) {
	    return 2;
	}

        //Use the factors to solve the linear system
        std::cout << "Solving the linear system...";
        Eigen::VectorXd w = solver.solve(r);

        // update for next iteration
        sol_vec_ = sol_vec_+w;
        for (int i=0; i<num_nodes_; i++) {
            f_vec.segment(i * k * dim_, k * dim_) =
                    f_vec.segment(i * k * dim_, k * dim_) + U.block(i * k * dim_, 0, k * dim_, dim_) *
                        w.block(i * dim_, 0, dim_, 1) +
                        p.segment(i * k * dim_, k * dim_);
        }

        if (w.norm() < ERR_TOL) return 0;
    }
    return 1;
}
