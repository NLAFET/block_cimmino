#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/extract.hpp>

#include "abcd.h"

using namespace boost::python;

BOOST_PYTHON_MODULE(abcdpy)
{

    class_<std::vector<int> >("_iVec")
        .def(vector_indexing_suite<std::vector<int> >()) ;
    class_<std::vector<double> >("_dVec")
        .def(vector_indexing_suite<std::vector<double> >()) ;

    class_<abcd>("abcd", init<>())
        .def("run", &abcd::operator())
        .def_readwrite("icntl", &abcd::icntl) 
        .def_readwrite("dcntl", &abcd::dcntl) 
        .def_readwrite("info", &abcd::info) 
        .def_readwrite("dinfo", &abcd::dinfo);
}
