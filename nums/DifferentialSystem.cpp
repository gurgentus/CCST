/* author: Gurgen Hayrapetyan
 * class representing a system of differential equations
 */

#include "DifferentialSystem.h"
#include <iostream>
#include <fstream>
#include <cassert>

/* Initialize with a pointer to the right-hand-side function, boundary condition function, and the
 * gradients of these functions.
 */
DifferentialSystem::DifferentialSystem(Eigen::VectorXd (*rhs)(double, const Eigen::VectorXd&),
                                       Eigen::MatrixXd (*rhs_grad)(double, const Eigen::VectorXd&),
                                       Eigen::VectorXd (*bc)(const Eigen::VectorXd &, const Eigen::VectorXd &),
                                       Eigen::MatrixXd (*bc_grad1)(const Eigen::VectorXd &, const Eigen::VectorXd &),
                                       Eigen::MatrixXd (*bc_grad2)(const Eigen::VectorXd &, const Eigen::VectorXd &),
                                       double t_min, double t_max, int dim)
{
    p_RhsFunc = rhs;
    p_RhsGradYFunc = rhs_grad;
    p_BcsFunc = bc;
    p_BcsGrad1Func = bc_grad1;
    p_BcsGrad2Func = bc_grad2;
    t_min_ = t_min;
    t_max_ = t_max;
    dim_ = dim;
}

/* Constructor method for classes inheriting this class.  The function pointers are automatically
 * pointed to corresponding class functions that will be overridden.
 */
DifferentialSystem::DifferentialSystem(double t_min, double t_max, int dim)
{
    t_min_ = t_min;
    t_max_ = t_max;
    dim_ = dim;
}

Eigen::VectorXd DifferentialSystem::RhsFunc(double t, const Eigen::VectorXd &y) {
    return ( *p_RhsFunc )(t, y);
}

Eigen::MatrixXd DifferentialSystem::RhsGradYFunc(double t, const Eigen::VectorXd &y) {
    return ( *p_RhsGradYFunc )(t, y);
}

Eigen::MatrixXd DifferentialSystem::BcsGrad1Func(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
    return ( *p_BcsGrad1Func )(y1, y2);
}

Eigen::MatrixXd DifferentialSystem::BcsGrad2Func(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
    return ( *p_BcsGrad2Func )(y1, y2);
}

Eigen::VectorXd DifferentialSystem::BcsFunc(const Eigen::VectorXd &y1, const Eigen::VectorXd &y2) {
    return ( *p_BcsFunc )(y1, y2);
}
