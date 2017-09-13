#ifndef NUMS_DIFFERENTIALSYSTEM_H
#define NUMS_DIFFERENTIALSYSTEM_H

#include "Eigen/Dense"

class DifferentialSystem {

private:
    // Eigen::MatrixXd coeffs_;     // Coefficient matrix of the ODE system
    Eigen::VectorXd (*p_RhsFunc)(double t, const Eigen::VectorXd& y);  // Function on RHS of ODE
    Eigen::MatrixXd (*p_RhsGradYFunc)(double t, const Eigen::VectorXd& y);  // Gradient with respect to y

    Eigen::VectorXd (*p_BcsFunc)(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Boundary condition function
    Eigen::MatrixXd (*p_BcsGrad1Func)(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Gradient wrt to first var
    Eigen::MatrixXd (*p_BcsGrad2Func)(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Gradient wrt to second var

protected:
    int dim_;   // State dimension

    // Interval for domain
    double t_min_;
    double t_max_;

public:
    double t_min() const {
        return t_min_;
    }

    double t_max() const {
        return t_max_;
    }

    int dim() const {
        return dim_;
    }

    DifferentialSystem(Eigen::VectorXd (*rhs)(double, const Eigen::VectorXd&),
                      Eigen::MatrixXd (*rhs_grad)(double, const Eigen::VectorXd&),
                      Eigen::VectorXd (*bc)(const Eigen::VectorXd& y1, const Eigen::VectorXd& y),
                      Eigen::MatrixXd (*bc_grad1)(const Eigen::VectorXd& y1, const Eigen::VectorXd& y),
                      Eigen::MatrixXd (*bc_grad2)(const Eigen::VectorXd& y1, const Eigen::VectorXd& y),
                      double t_min, double t_max, int dim);

    DifferentialSystem(double t_min, double t_max, int dim);

    virtual Eigen::VectorXd RhsFunc(double t, const Eigen::VectorXd& y);  // Function on RHS of ODE
    virtual Eigen::MatrixXd RhsGradYFunc(double t, const Eigen::VectorXd& y);  // Gradient with respect to y
    virtual Eigen::VectorXd BcsFunc(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Boundary condition function
    virtual Eigen::MatrixXd BcsGrad1Func(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Gradient wrt to first var
    virtual Eigen::MatrixXd BcsGrad2Func(const Eigen::VectorXd& left_bc, const Eigen::VectorXd& right_bc);  // Gradient wrt to second var

};

#endif //NUMS_DIFFERENTIALSYSTEM_H
