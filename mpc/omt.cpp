#include "Eigen/Dense"
#include "Eigen/StdVector"
#include <vector>
#include <iostream>
#include <boost/python/numpy.hpp>

using Eigen::MatrixXd;
using Eigen::VectorXd;

static const double ERR_TOL = 1e-8;
static const int MAX_ITER = 1000;

// orbital mechanics toolbox
class omt {

  std::vector<Eigen::Vector3d> r;
  std::vector<Eigen::Vector3d> v;

  Eigen::Vector3d r_res;
  Eigen::Vector3d v_res;

  double mu;
  double Omega;
  double e;
  double omega;
  double a;
  double period;
  double h;
  double i;

public:

  int lambert(boost::python::list& r1, boost::python::list& r2, double dt, bool prograde, const double mu)
  {
    r.resize(2);
    r[0] << boost::python::extract<double>(r1[0]), boost::python::extract<double>(r1[1]),
        boost::python::extract<double>(r1[2]);
    r[1] << boost::python::extract<double>(r2[0]), boost::python::extract<double>(r2[1]),
        boost::python::extract<double>(r2[2]);
    orbit_desc(r_res, v_res, r[0], r[1], dt, prograde, mu);
    // std::cout << r_res << std::endl;
    // std::cout << v_res << std::endl;
    return 1;
  }

  // boost::python::list get_v0()
  // {
  //   std::vector<double>::iterator iter;
	//   boost::python::list list;
  //   for (iter = v.begin(); iter != v.end(); ++iter) {
  //     list.append(*iter);
  //   }
  //   return list;
  // }

  boost::python::list get_v0()
  {
	  boost::python::list list;
    list.append(v_res(0));
    list.append(v_res(1));
    list.append(v_res(2));
    return list;
  }

  char const* greet()
  {
     return "hello, world";
  }

  /* Stumpff functions for use with universal variables */
  double stumpffS(double z)
  {
      double sq = sqrt(std::abs(z));
      if (z > 0) return (sq-sin(sq))/(sq*sq*sq);
      if (z < 0) return (sinh(sq)-sq)/(sq*sq*sq);
      return (1.0/6.0);
  }

  double stumpffC(double z)
  {
      if (z > 0) return (1-cos(sqrt(z)))/z;
      if (z < 0) return (cosh(sqrt(-z))-1)/(-z);
      return 0.5;
  }

  int orbit_desc(Eigen::Vector3d& r, Eigen::Vector3d& v, const Eigen::Vector3d& r1,
                      const Eigen::Vector3d& r2, double dt, bool prograde, const double mu)
  {
      double r1_norm = r1.norm();
      double r2_norm = r2.norm();
      Eigen::Vector3d C12 = r1.cross(r2); //QVector3D::crossProduct(r1, r2);
      double D12 = r1.dot(r2); //QVector3D::dotProduct(r1, r2);
      double P12 = r1_norm*r2_norm;
      double invcos = acos(D12/P12);
      double dtheta;
      if (prograde)
      {
          if (C12(2) >= 0)
          {
              dtheta = invcos;
          }
          else
          {
              dtheta = 2*M_PI-invcos;
          }
      }
      else
      {
          if (C12(2) < 0)
          {
              dtheta = invcos;
          }
          else
          {
              dtheta = 2*M_PI-invcos;
          }
      }
      double A = sin(dtheta)*sqrt(P12/(1-cos(dtheta)));

      double sqmu = sqrt(mu);

      double f, y, fprime, rat, pre, z, C, S;
      unsigned int iter;
      z = 0;
      for (iter = 0; iter < MAX_ITER; iter++) {

          C = stumpffC(z);
          S = stumpffS(z);
          pre = pow(y/C,1.5);
          y = r1_norm + r2_norm + A*(z*S-1)/sqrt(C);
          f = pre*S+A*sqrt(y)-sqmu*dt;
          if (fabs(z) > ERR_TOL)
          {
              fprime = pre*((1/(2*z))*(C-1.5*S/C)+0.75*S*S/C)+(A/8)*(3*S/C*sqrt(y)+A*sqrt(C/y));
          }
          else
          {
              fprime = (sqrt(2)/40)*pow(y,1.5)+(A/8)*(sqrt(y)+A*sqrt(0.5/y));
          }
          rat = f/fprime;
          if (std::abs(rat) < ERR_TOL) {
              break;
          }
          z = z - rat;
          // std::cout << f << ", " << fprime << ", " << rat << ", " << chi << std::endl;
      }
      if (iter == MAX_ITER) {
          return 1;
      }

      double lf = 1 - y/r1_norm;
      double lg = A*sqrt(y/mu);
      //double lfprime = (sqrt(mu)/P12)*sqrt(y/C)*(z*S-1);
      //double lgprime = 1 - y/r2_norm;
      Eigen::Vector3d v1 = (1/lg)*(r2-lf*r1);
      // QVector3D v2 = lgprime*r2/lg-(lf*lgprime-lfprime*lg)*r1/lg;

      v << v1(0), v1(1), v1(2);
      r << r1(0), r1(1), r1(2);

      orbit_desc(r, v, mu);
      return 0;
  }


  double get_mu() {
    return mu;
  }
  double get_Omega() {
    return Omega;
  }
  double get_e() {
    return e;
  }
  double get_omega() {
    return omega;
  }
  double get_a() {
    return a;
  }
  double get_period() {
    return period;
  }
  double get_h() {
    return h;
  }
  double get_i() {
    return i;
  }

  int orbit_desc(const Eigen::Vector3d& r, const Eigen::Vector3d& v, const double mu)
  {
    this->mu = mu;
    double r_scalar = r.norm();
    double speed = v.norm();
    double v_r = r.dot(v)/r_scalar;
    Eigen::Vector3d h_vec = r.cross(v);
    h = h_vec.norm();
    // inclination
    i = acos(h_vec.z()/h)*180/M_PI;
    // node line
    Eigen::Vector3d un;
    un << 0,0,1;
    Eigen::Vector3d N;
    N= un.cross(h_vec);
    double N_len = N.norm();
    // right ascension of the ascending node
    Omega = acos(N(0)/N_len);
    if (N(1) < 0) {
        Omega = 2*M_PI - Omega;
    }
    Omega = Omega*180/M_PI;
    // eccentricity vector
    Eigen::Vector3d e_vec = (1/mu)*((speed*speed - mu/r_scalar)*r - r_scalar*v_r*v);
    e = e_vec.norm();
    // argument of perigree
    omega = acos(N.dot(e_vec)/(N_len*e));
    if (e_vec(2) < 0) {
        omega = 2*M_PI - omega;
    }
    omega = omega*180/M_PI;

    secondary_params();
    return 0;
  }

  int secondary_params()
  {
      double r_p = h*h/(mu*(1+e));
      double r_a = h*h/(mu*(1-e));
      a = (r_p+r_a)/2;
      period = 2*M_PI*pow(a,1.5)/sqrt(mu);
      return 0;
  }

};

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(omt)
{
    // class_<std::vector<double> >("double_vector")
    //        .def(vector_indexing_suite<std::vector<double> >())
    //    ;
    //void (omt::*orbit_desc)(A&) = &Foo::orbit_desc;
    class_<omt>("omt")
    //  .def("show", &omt::show)
      .def("get_v0", &omt::get_v0)
      .def("greet", &omt::greet)
      .def("stumpffS", &omt::stumpffS)
      .def("stumpffC", &omt::stumpffC)
      .def("lambert", &omt::lambert)
      .def("get_mu", &omt::get_mu)
      .def("get_Omega", &omt::get_Omega)
      .def("get_e", &omt::get_e)
      .def("get_omega", &omt::get_omega)
      .def("get_a", &omt::get_a)
      .def("get_period", &omt::get_period)
      .def("get_h", &omt::get_h)
      .def("get_i", &omt::get_i)
    ;
}
