#include <boost/python.hpp>
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
#include <boost/python/numeric.hpp>
#include <boost/python/extract.hpp>
#include "boost/numpy.hpp"

#include "abcd.h"

namespace bp = boost::python;
namespace bn = boost::numpy;
using namespace std;

void set_matrix(
        abcd &obj,
        int const & m, int const & n, int const & nz,
        bn::ndarray const & irn,
        bn::ndarray const & jcn,
        bn::ndarray const & val
    ) 
{
    // check types
    if( (irn.get_dtype() != bn::dtype::get_builtin<int>()) ||
        (jcn.get_dtype() != bn::dtype::get_builtin<int>()) ||
        (val.get_dtype() != bn::dtype::get_builtin<double>())){
        PyErr_SetString(PyExc_TypeError, "Incorrect array data type");
        bp::throw_error_already_set();
    }
    
    obj.m = m;
    obj.n = n;
    obj.nz = nz;
    obj.irn = reinterpret_cast<int *>(irn.get_data());
    obj.jcn = reinterpret_cast<int *>(jcn.get_data());
    obj.val = reinterpret_cast<double *>(val.get_data());
}

void set_rhs(
    abcd &obj,
    bn::ndarray const & rhs,
    int const & nrhs)
{
    obj.nrhs = nrhs;
    obj.rhs = reinterpret_cast<double *>(rhs.get_data());
}

bn::ndarray get_sol(abcd &obj)
{
    return bn::from_data(obj.sol, bn::dtype::get_builtin<double>(),
                         bp::make_tuple(obj.n, obj.nrhs),
                         bp::make_tuple(sizeof(double), obj.n * sizeof(double)),
                         bp::object());
}


BOOST_PYTHON_MODULE(abcdpy)
{
    
    bn::initialize();  // have to put this in any module that uses Boost.NumPy

    bp::class_<std::vector<int> >("_iVec")
        .def(bp::vector_indexing_suite<std::vector<int> >()) ;
    bp::class_<std::vector<double> >("_dVec")
        .def(bp::vector_indexing_suite<std::vector<double> >()) ;

    bp::class_<abcd>("abcd", bp::init<>())
        .def("run", &abcd::operator())
        .def("set_matrix", set_matrix)
        .def("set_rhs", set_rhs)
        .def("get_sol", get_sol)
        .def_readwrite("icntl", &abcd::icntl) 
        .def_readwrite("dcntl", &abcd::dcntl) 
        .def_readwrite("info", &abcd::info) 
        .def_readwrite("dinfo", &abcd::dinfo);

    bp::enum_<Controls::icontrols>("icontrols")
        .value("nbparts", Controls::nbparts)
        .value("part_type", Controls::part_type)
        .value("part_guess", Controls::part_guess)
        .value("scaling", Controls::scaling)
        .value("itmax", Controls::itmax)
        .value("block_size", Controls::block_size)
        .value("verbose_level", Controls::verbose_level)
        .value("aug_type", Controls::aug_type)
        .value("aug_blocking", Controls::aug_blocking);
    
    bp::enum_<Controls::dcontrols>("dcontrols")
        .value("part_imbalance", Controls::part_imbalance)
        .value("threshold", Controls::threshold);

    bp::enum_<Controls::info>("info")
        .value("status", Controls::status)
        .value("nb_iter", Controls::nb_iter);
    
    bp::enum_<Controls::dinfo>("dinfo")
        .value("residual", Controls::residual)
        .value("forward_error", Controls::forward_error)
        .value("backward", Controls::backward)
        .value("scaled_residual", Controls::scaled_residual);
}
