#ifndef ABSTRACTDESOLVERDEF
#define ABSTRACTDESOLVERDEF

#include "DifferentialSystem.h"
#include "FiniteDifferenceGrid.h"
#include "Node.h"
#include <vector>
#include "Eigen/Dense"

class AbstractDeSolver
{
protected:
    int num_nodes_;  // number of grid nodes
    int dim_;  // dimension of the system
    FiniteDifferenceGrid grid_;  // pointer to a grid

    DifferentialSystem* p_ode_;  // pointer to a system to solve
    Eigen::VectorXd sol_vec_; // pointer to the solution vector
    Eigen::VectorXd f_vec; // y'
    std::string filename_;  // allow the user to specify the output file or use a default name

    static constexpr double ERR_TOL = 1e-1; // error tolerance for stopping iterations
    static const int MAX_ITER = 20; // maximum number of iterations before stopping

//protected:
//    double mStepSize;
//    double mInitialTime;
//    double mFinalTime;
//    double mInitialValue;
//    double t_ = 0;
//    std::vector<double> mInitialValueVector;
public:
    AbstractDeSolver(DifferentialSystem *p_ode, int num_nodes, Eigen::VectorXd& init_guess);
    virtual ~AbstractDeSolver() {};
    /* sets the filename if results should be written to a file */
    void set_filename(const std::string& name)
    {
        filename_ = name;
    }

    /* returns the solution vector */
    Eigen::VectorXd sol_vec() const
    {
        return sol_vec_;
    }

    /* returns grid point value */
    virtual double Step(int i) = 0;

    /* runs the solver */
    virtual int Solve() = 0;

    virtual void WriteSolutionFile() = 0;

//    void SetStepSize(double h);
//    void SetTimeInterval(double t0, double t1);
//    void SetInitialValue(double y0);
//    void SetInitialValue(std::vector<double> y0);
//    double time();
//
//    // virtual methods
//    virtual void InitialConditions() = 0;
//    virtual void UpdateState(double dt) = 0;
//    virtual void SolveEquation(std::vector<double> yi) = 0;

};

#endif
