#include <boost/python.hpp>
#include <boost/python/to_python_converter.hpp>
#include "qca.hpp"
#include "eigenHelpers.hpp"

namespace p = boost::python;

/*
 * Structs for automatic type conversion from C++ to Python.
 */
struct vector_to_python_list
{
    static PyObject* convert(std::vector<double> const& v)
    {
        p::list l;
        for (double d : v)
            l.append(d);
        return p::incref(l.ptr());
    }
};

struct vectorvector_to_python_list
{
    static PyObject* convert(std::vector<std::vector<double>> const& vv)
    {
        p::list k;
        for (const std::vector<double>& v : vv)
        {
            p::list l;
            for (double d : v)
                l.append(d);
            k.append(l);
        }
        return p::incref(k.ptr());
    }
};

struct dvector_to_python_list
{
    static PyObject* convert(DVector const& v)
    {
        p::list l;
        for (int i=0; i<v.size(); i++)
            l.append(v(i));
        return p::incref(l.ptr());
    }
};

/*
 * Template function to define a Python class for a QCA class.
 */
template<class QcaSystem>
void define_qca_python_class(const std::string& name)
{
    typedef double (QcaSystem::*Function1)(int) const;
    typedef std::vector<double> (QcaSystem::*Function2)(int) const;
    Function1 measurePolarization = &QcaSystem::measurePolarization;
    Function2 measureParticleNumber = &QcaSystem::measureParticleNumber;

    p::class_<QcaSystem>(name.c_str())
        .def(p::init<Layout>())
        .def("update",
             &QcaSystem::update,
             "Update the system for changed configuration.")
        .def("measurePolarization",
             measurePolarization,
             "Measure the polarization of a cell.")
        .def("measureParticleNumber",
             measureParticleNumber,
             "Measure the particle numbers of a cell.")
        .def("measureParticleNumberOverEnergy",
             &QcaSystem::measureParticleNumberOverEnergy,
             "Total particle number / occupancy per energy level.")
        .def("energies",
             &QcaSystem::energies,
             p::return_value_policy<p::copy_const_reference>(),
             "Retrieve the eigenenergies.")
        .def_readwrite("t",
                       &QcaSystem::t,
                       "Hopping parameter t.")
        .def_readwrite("V0",
                       &QcaSystem::V0,
                       "On-site Coulomb repulsion V_0.")
        .def_readwrite("beta",
                       &QcaSystem::beta,
                       "Inverse temperature beta.")
        .def_readwrite("l",
                       &QcaSystem::l,
                       "QCA cell layout.")
        ;
}


BOOST_PYTHON_MODULE (_qca)
{
    /*
     * This will enable user-defined docstrings and python signatures, while
     * disabling the C++ signatures
     */
    p::docstring_options my_docstring_options(true, true, false);

    /*
     * Register automatic type conversion from C++ to Python
     */
    p::to_python_converter<std::vector<double>, vector_to_python_list>();
    p::to_python_converter<std::vector<std::vector<double>>, vectorvector_to_python_list>();
    p::to_python_converter<const DVector, dvector_to_python_list>();

    /*
     * Layout
     */
    p::class_<Layout>("Layout", "QCA cell layout")
        .def("wire",
             &Layout::wire2e,
             "Construct a uniform wire.")
        .def("nonuniformWire",
             &Layout::nonuniformWire2e,
             "Construct a non-uniform wire.")
    ;

    /*
     * QCA systems
     */
    define_qca_python_class<QcaBond>("QcaBond");
    define_qca_python_class<QcaFixedCharge>("QcaFixedCharge");
    define_qca_python_class<QcaGrandCanonical>("QcaGrandCanonical");
}
