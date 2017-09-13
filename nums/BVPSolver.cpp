//
// Created by Gurgen Hayrapetyan on 9/9/17.
//

#include "BVPSolver.h"

#include <iostream>
#include <fstream>
#include <cassert>
#include "Eigen/Sparse"
#include "Eigen/SparseLU"

BVPSolver::BVPSolver(DifferentialSystem *p_ode, int num_nodes, int dim)
        : grid_(num_nodes, p_ode->t_min(), p_ode->t_max())
{
    p_ode_ = p_ode;
    dim_ = dim;
    num_nodes_ = num_nodes;
    sol_vec_ = Eigen::VectorXd::Zero(num_nodes_*dim_);
    filename_ = "ode_output.txt";
}

double BVPSolver::Step(int i)
{
    return grid_.nodes_[i].coordinate;
}

void BVPSolver::Solve()
{
    for (int iter=0; iter<MAX_ITER; iter++)
    {
        std::cout << "Iteration: " << iter << std::endl;
        Eigen::VectorXd u = Eigen::VectorXd(sol_vec_.segment(0,dim_));
        Eigen::VectorXd v = Eigen::VectorXd(sol_vec_.tail(dim_));

        Eigen::MatrixXd Ba = p_ode_->BcsGrad1Func(u,v);
        Eigen::MatrixXd Bb = p_ode_->BcsGrad2Func(u,v);
        Eigen::VectorXd beta = -p_ode_->BcsFunc(u,v);

        Eigen::SparseMatrix<double, Eigen::ColMajor> LeftMatrix(num_nodes_*dim_, num_nodes_*dim_);
        Eigen::VectorXd Q = Eigen::VectorXd::Zero(num_nodes_*dim_);
        for (int i=0; i<num_nodes_-1; i++)
        {
            double t = grid_.nodes_[i].coordinate;
            double tp = grid_.nodes_[i+1].coordinate;
            double h = tp-t;
            Eigen::VectorXd u = Eigen::VectorXd(sol_vec_.segment(dim_*i,dim_));
            Eigen::VectorXd up = Eigen::VectorXd(sol_vec_.segment(dim_*(i+1),dim_));

            Eigen::MatrixXd A = p_ode_->RhsGradYFunc(t, u);
            Eigen::MatrixXd Ap = p_ode_->RhsGradYFunc(tp, up);

            Eigen::MatrixXd I = Eigen::MatrixXd::Identity(dim_,dim_);
            Eigen::MatrixXd S = -(1.0/h)*I - 0.5*A;
            Eigen::MatrixXd R = (1.0/h)*I - 0.5*Ap;

            for (int j=0; j<dim_; j++)
            {
                for (int k=0; k<dim_; k++)
                {
                    LeftMatrix.coeffRef(dim_*i+j, dim_*i+k) += S(j,k);
                    LeftMatrix.coeffRef(dim_*i+j, dim_*(i+1)+k) += R(j,k);
                }
            }

            Q.segment(dim_*i,dim_) = 0.5*(p_ode_->RhsFunc(tp,up)+p_ode_->RhsFunc(t,u))-(up-u)/h;

        }

        for (int j=0; j<dim_; j++)
        {
            for (int k=0; k<dim_; k++)
            {
                LeftMatrix.coeffRef(dim_*(num_nodes_-1)+j, dim_+k) += Ba(j,k);
                LeftMatrix.coeffRef(dim_*(num_nodes_-1)+j, dim_*(num_nodes_-1)+k) += Bb(j,k);

            }
        }

        LeftMatrix.makeCompressed();
        Q.segment(dim_*(num_nodes_-1),dim_) = beta;

        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;

        // Compute the ordering permutation vector from the structural pattern of A
        solver.analyzePattern(LeftMatrix);
        // Compute the numerical factorization
        solver.factorize(LeftMatrix);
        //Use the factors to solve the linear system
        Eigen::VectorXd w = solver.solve(Q);

        sol_vec_ = sol_vec_+w;

        if (w.norm() < ERR_TOL) return;
    }
    return;
}

void BVPSolver::WriteSolutionFile()
{
    std::ofstream output_file(filename_.c_str());
    assert(output_file.is_open());
    for (int i=0; i<num_nodes_; i++)
    {
        double x = grid_.nodes_[i].coordinate;
        output_file << x << " " << sol_vec_(2*i) << "\n";
    }
    output_file.flush();
    output_file.close();
    std::cout << "Solution written to " << filename_ << std::endl;

}


//Eigen::VectorXd OrbitTransfer::RhsFunc(double t, const Eigen::VectorXd &y) {
//    Eigen::VectorXd rs = Eigen::VectorXd(2);
//    //rs << y(1), 34*sin(t)+4*y(0)-3*y(1);
//    rs << y(1), -9.8;
//    return rs;
//}
//
//Eigen::MatrixXd OrbitTransfer::RhsGradYFunc(double t, const Eigen::VectorXd &y) {
//    Eigen::MatrixXd rs = Eigen::MatrixXd(2,2);
////    rs << 0,1,
////            4,-3;
//    rs << 0, 1,
//            0, 0;
//    return rs;
//}
//
//Eigen::MatrixXd OrbitTransfer::BcsGrad1Func(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
//    Eigen::MatrixXd rs = Eigen::MatrixXd(2,2);
////    rs << 0,1,
////            0,0;
//    rs << 1,0,
//            0,0;
//    return rs;
//}
//
//Eigen::MatrixXd OrbitTransfer::BcsGrad2Func(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
//    Eigen::MatrixXd rs = Eigen::MatrixXd(2,2);
////    rs << 0,0,
////            1,0;
//    rs << 0,0,
//            1,0;
//    return rs;
//}
//
//Eigen::VectorXd OrbitTransfer::BcsFunc(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
//    Eigen::VectorXd rs = Eigen::VectorXd(2);
//    //rs << y1(1)+5, y2(0)-4;
//    rs << y1(0), y2(0)-40;
//    return rs;
//}