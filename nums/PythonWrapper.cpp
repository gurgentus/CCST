//
// Created by Gurgen Hayrapetyan on 9/12/17.
//

#include <boost/python.hpp>
#include <boost/python/numpy.hpp>
#include "OrbitTransfer.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(OrbitTransfer)
{
        // class_<std::vector<double> >("double_vector")
        //        .def(vector_indexing_suite<std::vector<double> >())
        //    ;
        //void (omt::*orbit_desc)(A&) = &Foo::orbit_desc;
        class_<OrbitTransfer>("OrbitTransfer")
                .def("greet", &OrbitTransfer::greet)
                .def("run", &OrbitTransfer::run)
        ;
}
