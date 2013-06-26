#include <boost/python.hpp>
#include <boost/python/to_python_converter.hpp>
#include <boost/python/enum.hpp>
#include "qca.hpp"
#include "eigenHelpers.hpp"

#include <memory>

namespace p = boost::python;

/*
 * Structs for automatic type conversion from C++ to Python.
 */
template<class T>
struct vector_to_python_list
{
    // static PyObject* convert(std::vector<double> const& v)
    // {
    //     p::list l;
    //     for (double d : v)
    //         l.append(d);
    //     return p::incref(l.ptr());
    // }
    
    static PyObject* convert(std::vector<T> const& v)
    {
        p::list l;
        for (T d : v)
            l.append(d);
        return p::incref(l.ptr());
    }
};

struct vector2d_to_python_tuple
{
    static PyObject* convert(Vector2d const& v)
    {
        p::tuple l = p::make_tuple(v(0), v(1));
        return p::incref(l.ptr());
    }
};

// struct vectorvector_to_python_list
// {
//     static PyObject* convert(std::vector<std::vector<double>> const& vv)
//     {
//         p::list k;
//         for (const std::vector<double>& v : vv)
//         {
//             p::list l;
//             for (double d : v)
//                 l.append(d);
//             k.append(l);
//         }
//         return p::incref(k.ptr());
//     }
// };

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
        .def("getLayout",
             &QcaSystem::getLayout,
             p::return_internal_reference<>(),
             "Retrieve the layout.")
        .def("setLayout",
             &QcaSystem::setLayout,
             "Set the layout.")
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
        .add_property("l_",
                      p::make_function(&QcaSystem::getLayout, p::return_internal_reference<>()),
                      &QcaSystem::setLayout)
        .def_readwrite("mu",
                       &QcaSystem::mu,
                       "Chemical potential mu.")
        ;
}

/*
 * Overloaded member functions in Layout
 */
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(wire_overloads, wire, 4, 5)
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(nonuniformWire_overloads, nonuniformWire, 4, 5)
// BOOST_PYTHON_MEMBER_FUNCTION_OVERLOADS(addDriverCell_overloads, addDriverCell, 4, 5)

void (Layout::*wire1)(int,double,double,double) = &Layout::wire;
void (Layout::*wire2)(int,double,double,double,ElectronsPerCell) = &Layout::wire;
void (Layout::*nonuniformWire1)(int,double,std::vector<double>,double) = &Layout::nonuniformWire;
void (Layout::*nonuniformWire2)(int,double,std::vector<double>,double,ElectronsPerCell) = &Layout::nonuniformWire;
void (Layout::*addDriverCell1)(double,double,double,double) = &Layout::addDriverCell;
void (Layout::*addDriverCell2)(double,double,double,double,ElectronsPerCell) = &Layout::addDriverCell;

struct LayoutWrap : Layout, p::wrapper<Layout>
{
    void addSite(double r_x, double r_y)
    {
        if (p::override addSite = this->get_override("addSite"))
            addSite(r_x,r_y);
        else
            Layout::addSite(r_x,r_y);
    }

    void default_addSite(double r_x, double r_y) { this->Layout::addSite(r_x,r_y); }
};

// class LayoutWrapper : public Layout
// {
// private:
//     PyObject* const self;
// public:
//     LayoutWrapper(PyObject* self_)
//     : self(self_)
//     {}
// 
//     // LayoutWrapper(PyObject* self_, const Layout& l)
//     // : Layout(l), self(self_)
//     // {}
// 
//     void addSite(double r_x, double r_y)
//     {
//         p::call_method<void>(self, "addSite", r_x, r_y);
//     }
// 
//     void default_addSite(double r_x, double r_y)
//     {
//         Layout::addSite(r_x,r_y);
//     }
// };

std::shared_ptr<Layout> test(const std::shared_ptr<Layout>& lw)
{
    //std::shared_ptr<Layout> l = std::make_shared<Layout>();
    std::shared_ptr<Layout> l; // = std::make_shared<Layout>(lw);
    l = lw;
    return l;
}

Layout* test2(Layout* lw)
{
    return lw;
}



/*
 * Define the qca module
 */
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
    //p::to_python_converter<std::vector<double>, vector_to_python_list>();
    p::to_python_converter<std::vector<double>, vector_to_python_list<double>>();
    p::to_python_converter<std::vector<Vector2d>, vector_to_python_list<Vector2d>>();
    p::to_python_converter<Vector2d, vector2d_to_python_tuple>();
    p::to_python_converter<std::vector<std::vector<double>>, vector_to_python_list<std::vector<double>>>();
    p::to_python_converter<const DVector, dvector_to_python_list>();

    /*
     * Layout
     */
    p::enum_<ElectronsPerCell>("ElectronsPerCell")
        .value("epc2", epc2)
        .value("epc6", epc6);
    //p::class_<Layout, std::shared_ptr<LayoutWrap>, boost::noncopyable>("Layout", "QCA cell layout")
    p::class_<LayoutWrap, boost::noncopyable>("Layout", "QCA cell layout")
    //p::class_<Layout, std::shared_ptr<Layout>, boost::noncopyable>("LayoutBase");
    //p::class_<Layout, LayoutWrapper, boost::noncopyable>("Layout", "QCA cell layout")
    //p::class_<Layout, std::shared_ptr<LayoutWrapper>, boost::noncopyable>("Layout", "QCA cell layout")
        // .def("wire",
        //      &Layout::wire,
        //      wire_overloads("Construct a uniform wire."))
        // .def("nonuniformWire",
        //      &Layout::nonuniformWire,
        //      nonuniformWire_overloads("Construct a non-uniform wire."))
        // .def("addDriverCell",
        //      &Layout::addDriverCell,
        //      addDriverCell_overloads("Add a driver cell.")[p::return_internal_reference<>()])
        .def("wire",
             wire1,
             "Construct a uniform wire.")
        .def("wire_",
             wire2,
             "Construct a uniform wire, with electrons-per-cell parameter.")
        .def("nonuniformWire",
             nonuniformWire1,
             "Construct a non-uniform wire.")
        .def("nonuniformWire_",
             nonuniformWire2,
             "Construct a non-uniform wire, with electrons-per-cell parameter.")
        .def("addDriverCell",
             addDriverCell1,
             "Add a driver cell.")
        .def("addDriverCell_",
             addDriverCell2,
             "Add a driver cell, with electrons-per-cell parameter.")
        // .def("addSite",
        //      &Layout::addSite,
        //      "Add a site.")
        .def("addSite",
             &Layout::addSite,
             &LayoutWrap::default_addSite,
             "Add a site.")
        .def("addCharge",
             &Layout::addCharge,
             "Add a charge.")
        .def("addCell",
             &Layout::addCell,
             "Add a cell.")
        .def("N_sites",
             &Layout::N_sites,
             "Number of sites.")
        .def("N_charges",
             &Layout::N_charges,
             "Number of charges.")
        .def("clear",
             &Layout::clear,
             "Clear / reset layout.")
        .add_property("r_sites", p::make_getter(&Layout::r_sites, p::return_value_policy<p::return_by_value>()))
        .add_property("r_charges", p::make_getter(&Layout::r_charges, p::return_value_policy<p::return_by_value>()))
        .add_property("charges", p::make_getter(&Layout::charges, p::return_value_policy<p::return_by_value>()))

    ;

    p::register_ptr_to_python<std::shared_ptr<Layout>>();
    //p::register_ptr_to_python<std::shared_ptr<LayoutWrapper>>();

    p::def("test", &test);
    p::def("test2", &test2, p::return_internal_reference<>());

    /*
     * QCA systems
     */
    define_qca_python_class<QcaBond>("QcaBond");
    define_qca_python_class<QcaFixedCharge>("QcaFixedCharge");
    define_qca_python_class<QcaGrandCanonical>("QcaGrandCanonical");
}
