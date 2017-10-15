/* author: Gurgen Hayrapetyan
 * car dynamics model and trajectory planners
 */

#include "Car.h"
#include <iostream>
#include <fstream>
#include "RungeKuttaSolver.h"
#include "Eigen/SparseLU"

Car::Car(double l, int N) : DifferentialSystem(0.0, 100.0, 5) {
    this->l = l;
    this->N = N;
    SetControlScheme(MinJerk);
}

Car::Car() : DifferentialSystem(0.0, 100.0, 5) {
    SetControlScheme(MinJerk);
    SetCarModel(SimpleCar);
}

void Car::SetParams(double l, double N)
{
    this->l = l;
    this->N = N;
}

void Car::SetCarModel(CarModel model)
{
    this->car_model = model;
}

void Car::SetControlScheme(ControlScheme scheme)
{
    this->control_scheme = scheme;
}

/* Dynamics of the Car */
Eigen::VectorXd Car::RhsFunc(double t, const Eigen::VectorXd& y_vec)
{
    Eigen::VectorXd rs = Eigen::VectorXd(5);
    double x = y_vec(0);
    double y = y_vec(1);
    double v = y_vec(2);
    double xi = y_vec(3);
    double m = y_vec(4);

    double D = 0;
    double g = 0;

    UpdateControls(y_vec, MinTime);

    if (car_model == SimpleCar) {
        rs <<   v*cos(xi),
                v*sin(xi),
                (T_*cos(dw_)*cos(dw_)-D)/m,
                (v/l)*tan(dw_),
                g;
    }
    return rs;
}

/* Sets the control variables based on the control scheme */
void Car::UpdateControls(const Eigen::VectorXd& y_vec, ControlScheme control_scheme)
{
    if (control_scheme == MinTime) {
        dw_ = 0;
        T_ = MinTimeStoppageThrottleControl(y_vec);
    }
}

/* Minimum time stoppage control scheme based on the maximum principle derived
 * bang-bang control scheme */
double Car::MinTimeStoppageThrottleControl(const Eigen::VectorXd& y_vec)
{
    double x = y_vec(0);
    double y = y_vec(1);
    double v = y_vec(2);
    double xi = y_vec(3);
    double m = y_vec(4);

    double dw = 0;
    double T = -1;
    if (x < -0.5*fabs(v)*v) {
        T = 1;
    }
    if (fabs(v) < 0.001) {
        T = 0;
    }
    return T;
}

/* Gradient of the dynamics equation for faster calculations in some numerical schemes */
Eigen::MatrixXd Car::RhsGradYFunc(double t, const Eigen::VectorXd& y_vec)
{
    double x = y_vec(0);
    double y = y_vec(1);
    double v = y_vec(2);
    double xi = y_vec(3);
    double m = y_vec(4);

    Eigen::MatrixXd rs = Eigen::MatrixXd(5,5);
    rs << 0,0,cos(xi),-v*sin(xi),0,
            0,0,sin(xi),v*cos(xi),0,
            0,0,0,0,0,
            0,0,tan(dw_)/l,0,0,
            0,0,0,0,0;

    return rs;
}

/* Boundary conditions */
Eigen::VectorXd Car::BcsFunc(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
    double x = y1(0);
    double y = y1(1);
    double v = y1(2);
    double xi = y1(3);
    double m = y1(4);

    Eigen::VectorXd rs = Eigen::VectorXd(5);
    rs << x-this->x, y-this->y, v-this->v, xi-this->xi, m-this->m;
    return rs;
}

/* Gradients of the boundary conditions for faster calculations in some numerical schemes */
Eigen::MatrixXd Car::BcsGrad1Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
    Eigen::MatrixXd rs = Eigen::MatrixXd(5,5);
    rs <<   1,0,0,0,0,
            0,1,0,0,0,
            0,0,1,0,0,
            0,0,0,1,0,
            0,0,0,0,1;

    return rs;
}

Eigen::MatrixXd Car::BcsGrad2Func(const Eigen::VectorXd& y1, const Eigen::VectorXd& y2)
{
    Eigen::MatrixXd rs = Eigen::MatrixXd(5,5);
    rs <<   0,0,0,0,0,
            0,0,0,0,0,
            0,0,0,0,0,
            0,0,0,0,0,
            0,0,0,0,0;

    return rs;
}

/* Trajectory generator based on minimizing the square of the given derivative of the position between given waypoints.
 * For example, setting num_coeff = 6, and only two waypoints (start and end position) generates a jerk minimizing
 * trajectory between the two points (note enough conditions have to be supplied at the waypoints).
 * As long as the number of conditions is consistent they can be specified more generally, to get more general behavior
 * than using a spline.
 */
int Car::GenerateMinDerivTraj(int num_coeffs, std::vector<Eigen::MatrixXd>& waypoints, Eigen::VectorXd& time_at_waypoint, Eigen::MatrixXd& spline_coeffs)
{
    int dims = waypoints.size();
    int n = (int)waypoints[0].cols();
    // define spline coefficients based on degree, e.g. num_coeffs = 6 is jerk minimizing
    // num_coeffs = 8 is snap minimizing etc.
    spline_coeffs = Eigen::MatrixXd::Zero(dims,num_coeffs*(n-1));
    int conds = num_coeffs/2;

    // Variational problem solutions are polynomyals of degree num_coeffs-1
    // We setup a linear system corresponding to the boundary and possibly
    // midpoint conditions given as part of waypoints matrix
    for (int dim = 0; dim<dims; dim++) {
        // This will be a large sparse matrix so use Eigen::SparseMatrix
        Eigen::SparseMatrix<double, Eigen::ColMajor> A(num_coeffs*(n-1), num_coeffs*(n-1));
        Eigen::VectorXd b = Eigen::VectorXd::Zero(num_coeffs*(n-1));

        for (int i = 0; i < n - 1; i++) { // condition at T_i, where T_0 = 0
            for (int j = 0; j < conds; j++) { // condition on the j'th derivative
                A.coeffRef(conds * i + j, num_coeffs * i + j) = Factorial(j);
                if (i == 0) {
                    b(j) = waypoints[dim](j,0);
                }
            }
        }

        for (int i = 1; i < n; i++) { // condition at T_i, where T_0 = 0
            double time_interval = time_at_waypoint(i) - time_at_waypoint(i - 1);
            double prefix = -1;
            if (i == n-1)
            {
                prefix = 1;
            }
            for (int j = 0; j < conds; j++) { // condition on the j'th derivative
                // j'th derivative at time T_i, this is conds*i+j row of the matrix
                for (int k = j; k < num_coeffs; k++) {
                    A.coeffRef(conds * i + j, num_coeffs * (i - 1) + k) = A.coeffRef(conds * i + j, num_coeffs * (i - 1) + k) +
                         prefix * (Factorial(k) / Factorial(k - j)) * pow(time_interval, k - j);
                }
                if (i == n - 1) {
                    b(conds*i + j) = waypoints[dim](j,n-1);
                }
            }
        }
        for (int i = 1; i < n-1; i++) {
            double time_interval = time_at_waypoint(i) - time_at_waypoint(i - 1);
            // at each intermediate value need conds-1 more derivative conditions and waypoint value condition
            for (int j = 0; j < conds-1; j++) { // condition on the j'th derivative
                A.coeffRef(conds * n + j, num_coeffs * i + j) = Factorial(j);
                // j'th derivative at time T_i, this is conds*i+j row of the matrix
                for (int k = conds+j; k < num_coeffs; k++) {
                    A.coeffRef(conds * n + j, num_coeffs * (i - 1) + k) =
                            A.coeffRef(conds * n + j, num_coeffs * (i - 1) + k) -
                                    (Factorial(k) / Factorial(k - conds- j)) * pow(time_interval, k - conds - j);
                }
            }
            A.coeffRef(conds * (n + 1) - 2 + i, num_coeffs * i) = 1;
            b(conds * (n + 1) - 2 + i) = waypoints[dim](0, i); // position constraints
        }


        // prepare for solving the system with sparse matrix
        A.makeCompressed();
        Eigen::SparseLU<Eigen::SparseMatrix<double, Eigen::ColMajor>> solver;
        // Compute the ordering permutation vector from the structural pattern of Gamma
        solver.analyzePattern(A);
        // Compute the numerical factorization
        solver.factorize(A);

        std::cout << "mat: " << Eigen::MatrixXd(A) << std::endl;
        std::cout << "b: " << b << std::endl;
        std::cout << "Solving the linear system...";

        std::cout << "Solver info: " << solver.info() << std::endl;
        if (solver.info() == 1) {
            return 2;
        }
        //Use the factors to solve the linear system
        spline_coeffs.row(dim) = solver.solve(b).transpose();

    }
    std::cout << spline_coeffs << std::endl;

    return 0;


}

/* time spent between waypoints defaults to proportional to distance between waypoints */
int Car::GenerateMinDerivTraj(const Eigen::MatrixXd& waypoints, Eigen::VectorXd& spline_coeffs)
{
    // TODO: setup time_between based on distance and call GenerateMinDerivTraj:
    return 0;
}

/* Calculate position along the trajectory based on trajectory polynomial coefficients */
void Car::DesiredState(int num_coeffs, double t, Eigen::MatrixXd& spline_coeffs, const Eigen::VectorXd& time_at_waypoint, Eigen::VectorXd& state) {
    int t_index = 0;
    if (t > time_at_waypoint(time_at_waypoint.size() - 1)) {
        t_index = (int) time_at_waypoint.size() - 1;
    }
    for (int i = 0; i < time_at_waypoint.size() - 1; i++) {
        if (t < time_at_waypoint(i + 1) && t > time_at_waypoint(i)) {
            t_index = i;
        }
    }

    Eigen::MatrixXd coeff = spline_coeffs.block(0, t_index * num_coeffs, state.size(), num_coeffs);
    state = Eigen::VectorXd::Zero(state.size());
    if (t == 0)
    {
        state = coeff.col(0);
    }
    else
    {
        for (int i = 0; i < num_coeffs; i++) {
            state = state + pow((t - time_at_waypoint(t_index)),i) * coeff.col(i);
        }
    }
}

/* Test code for simulating min jerk lane change */
int Car::SimulateMinJerkLaneChange(const Eigen::VectorXd& init_y, std::vector<Eigen::MatrixXd>& waypoints, Eigen::VectorXd& time_at_waypoint)
{
    std::cout << "Starting solution simulation" << std::endl;

    Eigen::VectorXd init_state = Eigen::VectorXd::Zero(dim());
    x = init_y(0);
    y = init_y(1);
    v = init_y(2);
    xi = init_y(3);
    m = init_y(4);


    Eigen::MatrixXd spline_coeffs;

    int status = GenerateMinDerivTraj(6, waypoints, time_at_waypoint, spline_coeffs);
    if (status == 2) {
        return 2;
    }

    std::ofstream output_file("traj_history.txt");
    assert(output_file.is_open());
    double t = 0.0;
    Eigen::VectorXd state = Eigen::VectorXd::Zero(2);
    for (int i=0; i<N; i++)
    {
        DesiredState(6, t, spline_coeffs, time_at_waypoint, state);

        output_file << state(0) << " " << state(1) << "\n";
        t = t + time_at_waypoint(time_at_waypoint.size()-1)/N;
    }
    output_file.flush();
    output_file.close();
    std::cout << "Solution written to traj_history.txt" << std::endl;

    return 0;
}

/* Test code for simulating min time stopping */
int Car::SimulateMinTimeBreaking(const Eigen::VectorXd& init_y)
{
    std::cout << "Starting solution simulation" << std::endl;

    Eigen::VectorXd init_state = Eigen::VectorXd::Zero(dim());
    x = init_y(0);
    y = init_y(1);
    v = init_y(2);
    xi = init_y(3);
    m = init_y(4);

    Eigen::VectorXd init_guess = Eigen::VectorXd::Zero(dim());
    init_guess(0) = x;
    init_guess(1) = y;
    init_guess(2) = v;
    init_guess(4) = m;

    // Instantiate a CollocationSolver with 1000 grid points for a system of
    RungeKuttaSolver ode(this, N, init_guess);
    std::cout << "Solving..." << std::endl;

    int status = ode.Solve();
    if (status == 0) {
        std::cout << "Solved." << std::endl;
        sol_vec_ = ode.sol_vec();
        std::cout << sol_vec_.tail(dim()) << std::endl;
    }

    std::cout << "Solved." << std::endl;
    ode.WriteSolutionFile();

    std::ofstream output_file("control_history.txt");
    assert(output_file.is_open());
    double t = 0.0;
    for (int i=0; i<N; i++)
    {
        output_file << t << " " << MinTimeStoppageThrottleControl(sol_vec_.segment(i*dim(),dim())) << "\n";
        t = t + 0.01;
    }
    output_file.flush();
    output_file.close();
    std::cout << "Solution written to control_history.txt" << std::endl;

    return 0;
}

/* utility function */
int Car::Factorial(int n)
{
    return (n == 1 || n == 0) ? 1 : Factorial(n - 1) * n;
}

//boost::python::list OrbitTransfer::getT()
//{
//    boost::python::list result;
//    double h = 1.0/N;
//    for (int i=0; i<N; i++)
//    {
//        result.append(i*h*tf);
//    }
//    return result;
//}
//
//boost::python::list OrbitTransfer::getAngle()
//{
//    boost::python::list result;
//    for (int i=0; i<N; i++)
//    {
//        result.append(180*atan(sol_vec_(i*dim_+5)/sol_vec_(i*dim_+6))/3.14);
//    }
//    return result;
//}
//

boost::python::list Car::GenerateOptimalStoppingTrajectory(boost::python::list& start_pos, boost::python::list& start_vel)
{
    boost::python::list result;

    double cur_posx = boost::python::extract<double>(start_pos[0]);
    double cur_posy = boost::python::extract<double>(start_pos[1]);
    double cur_velx = boost::python::extract<double>(start_vel[0]);
    double cur_vely = boost::python::extract<double>(start_vel[1]);

    double car_speed = sqrt(cur_velx*cur_velx + cur_vely*cur_vely);

    Eigen::VectorXd init_state = Eigen::VectorXd::Zero(5);
    init_state << cur_posx, cur_posy, car_speed, 0, 1;


    // Instantiate a CollocationSolver with 1000 grid points for a system of
    RungeKuttaSolver ode(this, N, init_state);
    std::cout << "Solving..." << std::endl;

    int status = ode.Solve();
    if (status == 0) {
        std::cout << "Solved." << std::endl;
        sol_vec_ = ode.sol_vec();
        std::cout << sol_vec_.tail(dim()) << std::endl;
    }

    std::cout << "Solved." << std::endl;

    Eigen::VectorXd state = Eigen::VectorXd::Zero(2);
    double t = 0;
    for (int i=0; i<N; i++)
    {
        result.append(t);
        result.append(MinTimeStoppageThrottleControl(sol_vec_.segment(i*dim(),dim())));
        t = t + 0.01;
    }

    return result;
}

boost::python::list Car::Generate2DTrajectory(boost::python::list& start_pos, boost::python::list& start_vel,
                                              boost::python::list& start_acc, boost::python::list& end_pos,
                                              boost::python::list& end_vel, boost::python::list& end_acc, double time_ahead, bool fuzzy, int num)
{
    double cur_posx = boost::python::extract<double>(start_pos[0]);
    double cur_posy = boost::python::extract<double>(start_pos[1]);
    double cur_velx = boost::python::extract<double>(start_vel[0]);
    double cur_vely = boost::python::extract<double>(start_vel[1]);
    double cur_accx = boost::python::extract<double>(start_acc[0]);
    double cur_accy = boost::python::extract<double>(start_acc[1]);
    double car_speed = sqrt(cur_velx*cur_velx + cur_vely*cur_vely);

    double ref_posx = boost::python::extract<double>(end_pos[0]);
    double ref_posy = boost::python::extract<double>(end_pos[1]);
    double ref_velx = boost::python::extract<double>(end_vel[0]);
    double ref_vely = boost::python::extract<double>(end_vel[1]);
    double ref_accx = boost::python::extract<double>(end_acc[0]);
    double ref_accy = boost::python::extract<double>(end_acc[1]);

    boost::python::list results;

    std::normal_distribution<double> ref_posx_distr(ref_posx,0.3);
    std::normal_distribution<double> ref_velx_distr(ref_velx,0.3);
    std::normal_distribution<double> ref_posy_distr(ref_posy,0.3);
    std::normal_distribution<double> ref_vely_distr(ref_vely,0.3);

    double h = time_ahead/N;
    int num_feasible_traj = 0;
    double best_cost = 10000;

    for (int i=0; i<num; i++) {
        boost::python::list x_result;
        boost::python::list y_result;
        std::vector<double> x_vec;
        std::vector<double> y_vec;

        if (fuzzy) {
            ref_posx = ref_posx_distr(generator);
            ref_vely = ref_velx_distr(generator);
            ref_posy = ref_posy_distr(generator);
            ref_vely = ref_vely_distr(generator);
        }

        Eigen::VectorXd init_state = Eigen::VectorXd::Zero(5);
        init_state << cur_posx, cur_posy, car_speed, 0, 1;

        std::vector<Eigen::MatrixXd> waypoints(2);
        waypoints[0] = Eigen::MatrixXd(3,2);
        waypoints[0] << cur_posx, ref_posx,
                cur_velx, ref_velx,
                cur_accx, ref_accx;
        waypoints[1] = Eigen::MatrixXd(3,2);
        waypoints[1] << cur_posy, ref_posy,
                cur_vely, ref_vely,
                cur_accy, ref_accy;
        std::cout << waypoints[0] << std::endl << waypoints[1] << std::endl;
        Eigen::VectorXd time_at_waypoint = Eigen::VectorXd(2);
        time_at_waypoint << 0, time_ahead;

        Eigen::MatrixXd spline_coeffs;

        int status = GenerateMinDerivTraj(6, waypoints, time_at_waypoint, spline_coeffs);
        if (status == 2) {
            std::cout << "error" << std::endl;
        }
        Eigen::VectorXd state = Eigen::VectorXd::Zero(2);
        double t = 0;
        bool feasible_trajectory = true;
        double total_cost = 0;

        for (int i=0; i<N; i++)
        {
            DesiredState(6, t, spline_coeffs, time_at_waypoint, state);
            x_result.append(state(0));
            y_result.append(state(1));
            x_vec.push_back(state(0));
            y_vec.push_back(state(1));
            t = t + h;
            double x_cost = RunningCost(x_vec, h, true, true);
            if (x_cost > 10.0) {
                feasible_trajectory = false;
                break;
            }
            double y_cost = RunningCost(y_vec, h, true, true);
            if (y_cost > 5.0) {
                feasible_trajectory = false;
                break;
            }
            total_cost = total_cost + x_cost + y_cost;
        }
        if (feasible_trajectory) {
            num_feasible_traj++;
            if (total_cost < best_cost) {
                results.insert(0,x_result);
                results.insert(1,y_result);
                best_cost = total_cost;
            }
            else {
                results.append(x_result);
                results.append(y_result);
            }
        }

    }
    return results;
}

double Car::RunningCost(std::vector<double>& state, double h, bool dir_check, bool jerk_check)
{
    double cost = 0;
    int n = state.size();
    if (dir_check) {
        if (n > 2) {
            if (state[n-1] < state[n-2]) {
                cost = cost + 100.0;
            }
        }
    }
    if (jerk_check) {
        if (n > 3) {
            cost = cost + (state[n-1]-3*state[n-2]+3*state[n-3]-state[n-4])/(h*h*h);
        }
    }
    return cost;
}

#include <boost/python.hpp>
//#include <boost/python/numpy.hpp>

using namespace boost::python;

BOOST_PYTHON_MODULE(Car)
{
    class_<Car>("Car")
//            .def("greet", &Car::greet)
            .def("GenerateOptimalStoppingTrajectory", &Car::GenerateOptimalStoppingTrajectory)
            .def("Generate2DTrajectory", &Car::Generate2DTrajectory)
            .def("SetParams", &Car::SetParams)
            ;
}
