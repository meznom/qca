#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cqca test
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/json_parser.hpp>

#define STORAGE_TYPE_OF_FERMIONIC_STATE uint32_t
#include "cqca.hpp"

using boost::property_tree::ptree;
using namespace boost::property_tree::json_parser;

std::string jsonify (std::string s)
{
    size_t pos = 0;
    while ((pos=s.find("'")) != std::string::npos)
    {
        s.replace(pos, 1, "\"");
    }
    return s;
}

ptree ptreeFromJson (const std::string& s)
{
    ptree p;
    std::stringstream ss(jsonify(s));
    read_json(ss, p);
    return p;
}

std::string ptreeToJson (const ptree& p)
{
    std::stringstream ss;
    write_json(ss, p);
    return ss.str();
}

template<typename T>
std::string toString (T v)
{
    std::stringstream ss;
    ss << std::setprecision(16) << v;
    return ss.str();
}

bool epsilonEqual (double v, double w, double epsilon = 10E-10)
{
    return fabs(v-w) < epsilon;
}

BOOST_AUTO_TEST_CASE ( test_configurable_layout )
{
    CLayout l1;
    l1.setConfig(ptreeFromJson("{}")); // default configuration
    BOOST_CHECK (l1.getConfig() == ptreeFromJson(
                "{'type': 'wire', 'cells': 1, 'a': 1, 'b': 3, 'Pext': 0, 'epc': 2}"));

    CLayout l2;
    BOOST_CHECK_THROW (
        l2.setConfig(ptreeFromJson(
                "{'type': 'wireblah', 'cells': 1, 'a': 1, 'b': 3, 'Pext': 0, 'epc': 2}")),
        ConfigurationException
        );

    CLayout l3;
    BOOST_CHECK_THROW (
        l3.setConfig(ptreeFromJson(
                "{'type': 'wire', 'cells': -2, 'a': 1, 'b': 3, 'Pext': 0, 'epc': 2}")),
        ConfigurationException
        );

    CLayout l4;
    BOOST_CHECK_THROW (
        l4.setConfig(ptreeFromJson(
                "{'type': 'wire', 'cells': 3, 'a': 1, 'b': -3.1, 'Pext': 0, 'epc': 2}")),
        ConfigurationException
        );

    CLayout l5;
    BOOST_CHECK_THROW (
        l5.setConfig(ptreeFromJson(
                "{'type': 'wire', 'cells': 3, 'a': 1, 'b': 3.1, 'Pext': 2.1, 'epc': 2}")),
        ConfigurationException
        );

    CLayout l6;
    BOOST_CHECK_THROW (
        l6.setConfig(ptreeFromJson(
                "{'type': 'wire', 'cells': 3, 'a': 1, 'b': 3.1, 'Pext': 0.1, 'epc': 4}")),
        ConfigurationException
        );

    CLayout l8;
    l8.setConfig(ptreeFromJson(
                "{'type': 'wire', 'cells': 4, 'a': 2.1, 'b': 3.5, 'Pext': 0.3, 'epc': 6}"));
    BOOST_CHECK (l8.getConfig() == ptreeFromJson(
                "{'type': 'wire', 'cells': 4, 'a': 2.1, 'b': 3.5, 'Pext': 0.3, 'epc': 6}"));
    BOOST_CHECK (l8.layout().N_sites() == 16);
    BOOST_CHECK (l8.layout().N_charges() == 4);

    CLayout l9;
    l9.setConfig(ptreeFromJson(
                "{'type': 'nonuniformwire', 'cells': 3, 'a': 0.1, 'bs': [1,1,1.2], 'Pext': 0, 'epc': 2}"));
    BOOST_CHECK (l9.getConfig() == ptreeFromJson(
                "{'type': 'nonuniformwire', 'cells': 3, 'a': 0.1, 'bs': [1,1,1.2], 'Pext': 0, 'epc': 2}"));
    Layout l9_;
    std::vector<double> bs;
    bs.push_back(1);
    bs.push_back(1);
    bs.push_back(1.2);
    l9_.addNonuniformWire(0,0, 3, 0.1, bs, 0);
    BOOST_CHECK (l9.layout() == l9_);
}

BOOST_AUTO_TEST_CASE ( test_configurable_qca_systems )
{
    CQca s1;
    BOOST_CHECK (s1.getConfig() == ptree());
    s1.measure();
    // operator== does not work for rtree
    //BOOST_CHECK (s1.measure() == rtree());

    CQca s2;
    s2.setConfig(ptreeFromJson("{}")); // default configuration
    BOOST_CHECK (s2.getConfig() == ptreeFromJson(
        "{'model': 'grandcanonical', 't': 1, 'td': 0, 'Vext': 0, 'V0': 1000, "
        "'mu': 0, 'epsilonr': " + toString(QCA_NATURAL_EPSILON_R) + 
        ", 'lambdaD': 0, 'q': 0, 'beta': 1, "
        "'layout': {'type': 'wire', 'cells': 1, 'a': 1, "
        "'b': 3, 'Pext': 0, 'epc': 2}, 'observables': {}}"));

    CQca s3;
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 3}, 'observables': {'P': 'all2'}}"));
    BOOST_CHECK_THROW(s3.measure(), ConfigurationException);
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 3}, 'observables': {'P': 4}}"));
    BOOST_CHECK_THROW(s3.measure(), ConfigurationException);
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 2}, 'observables': {'P': [1,2]}}"));
    BOOST_CHECK_THROW(s3.measure(), ConfigurationException);
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 2}, 'observables': {'P': [1,-1]}}"));
    BOOST_CHECK_THROW(s3.measure(), ConfigurationException);
    // these should not throw anything
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 2}, 'observables': {'P': 'all'}}"));
    s3.measure();
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 2}, 'observables': {'P': [1,0]}}"));
    s3.measure();
    s3.setConfig(ptreeFromJson(
        "{'model': 'bond', 'layout': {'cells': 3}, 'observables': {'P': [1,2,0]}}"));
    s3.measure();

    CQca s4;
    s4.setConfig(ptreeFromJson(
        "{'model': 'bond', 't': 2, 'td': 0.1, 'Vext': 0.001, 'V0': 500, "
        "'mu': -2, 'epsilonr': 1.3, 'lambdaD': 6.5, 'q': 1, 'beta': 10, "
        "'layout': {'type': 'nonuniformwire', 'cells': 3, 'a': 1.2, "
        "'bs': [1.2,1,2], 'Pext': 0.7, 'epc': 6}, 'observables': {'P': 'all'}}"));
    BOOST_CHECK (s4.getConfig() == ptreeFromJson(
        "{'model': 'bond', 't': 2, 'td': 0.1, 'Vext': 0.001, 'V0': 500, "
        "'mu': -2, 'epsilonr': 1.3, 'lambdaD': 6.5, 'q': 1, 'beta': 10, "
        "'layout': {'type': 'nonuniformwire', 'cells': 3, 'a': 1.2, "
        "'bs': [1.2,1,2], 'Pext': 0.7, 'epc': 6}, 'observables': {'P': 'all'}}"));
    
    double a = 1.0/250.0;
    double b = 1.75 * a;
    double V0 = 10 / a;
    s4.setConfig(ptreeFromJson(
        "{'model': 'bond', 't': 1, 'td': 0.2, 'V0': " + toString(V0) + ", 'beta': 1000000, "
        "'layout': {'type': 'wire', 'cells': 1, 'a': " + toString(a) + ", "
        "'b': " + toString(b) + ", 'Pext': 0.01}, 'observables': {'P': 'all'}}"));
    BOOST_CHECK (s4.getConfig() == ptreeFromJson(
        "{'model': 'bond', 't': 1, 'td': 0.2, 'Vext': 0, 'V0': " + toString(V0) + ", "
        "'mu': 0, 'epsilonr': " + toString(QCA_NATURAL_EPSILON_R) + ", 'lambdaD': 0, "
        "'q': 0, 'beta': 1000000, "
        "'layout': {'type': 'wire', 'cells': 1, 'a': " + toString(a) + ", "
        "'b': " + toString(b) + ", 'Pext': 0.01, 'epc': 2}, 'observables': {'P': 'all'}}"));
    rtree r = s4.measure();
    BOOST_CHECK (epsilonEqual(r.get<double>("P.0"), 0.29, 0.01));

    s4.setConfig(ptreeFromJson(
        "{'model': 'bond', 't': 1, 'td': 0, 'V0': 1000, 'beta': 1, "
        "'layout': {'type': 'wire', 'cells': 3, 'a': 0.01, "
        "'b': 0.023, 'Pext': 1}, 'observables': {'P': 'all'}}"));
    r = s4.measure();
    BOOST_CHECK (epsilonEqual(r.get<double>("P.0"), 0.670243, 1E-5));
    BOOST_CHECK (epsilonEqual(r.get<double>("P.1"), 0.462795, 1E-5));
    BOOST_CHECK (epsilonEqual(r.get<double>("P.2"), 0.295182, 1E-5));
    
    a = 1.0/250.0;
    b = 5 * a;
    V0 = 10.0 / a;
    s4.setConfig(ptreeFromJson(
        "{'model': 'fixed', 't': 1, 'td': 0.2, 'V0': " + toString(V0) + ", " 
        "'beta': 1000, 'layout': {'type': 'wire', 'cells': 2, "
        "'a': " + toString(a) + ", 'b': " + toString(b) + ", 'Pext': 0.1, "
        "'epc': 6}, 'observables': {'P': 'all', 'N': 0}}"));
    r = s4.measure();
    BOOST_CHECK (epsilonEqual(r.get<double>("N.0.total"), 6));
    BOOST_CHECK (epsilonEqual(r.get<double>("P.0"), 0.268, 0.001));
    BOOST_CHECK (epsilonEqual(r.get<double>("P.1"), 0.212, 0.001));
    
    //std::cerr << ptreeToJson(s3.getConfig()) << std::endl;
}

BOOST_AUTO_TEST_CASE ( test_vparam )
{
    ptree p;
    p = ptreeFromJson("{'a': 'blah'}");
    BOOST_CHECK_THROW (
        VParam v1("v1", p.get_child("a")),
        ConfigurationException);
    
    p = ptreeFromJson("{'a': '3,4,5'}");
    BOOST_CHECK_THROW (
        VParam v2("v2", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': '3;4'}");
    BOOST_CHECK_THROW (
        VParam v3("v3", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': '3;;'}");
    BOOST_CHECK_THROW (
        VParam v4("v4", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': '1;2;-0.2'}");
    BOOST_CHECK_THROW (
        VParam v5("v5", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': '2.1;1.5;0.1'}");
    BOOST_CHECK_THROW (
        VParam v6("v6", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': '1;2;0.b5'}");
    BOOST_CHECK_THROW (
        VParam v7("v7", p.get_child("a")),
        ConfigurationException);
    
    p = ptreeFromJson("{'a': '1;2a;0.1'}");
    BOOST_CHECK_THROW (
        VParam v8("v8", p.get_child("a")),
        ConfigurationException);

    p = ptreeFromJson("{'a': [1,2,3,4,5]}");
    VParam v9("v9", p.get_child("a"));
    BOOST_CHECK (v9.size() == 5);
    BOOST_CHECK (v9.index == 0);
    BOOST_CHECK (v9.currentValue() == "1");
    v9.index++;
    BOOST_CHECK (v9.currentValue() == "2");
    
    p = ptreeFromJson("{'a': '1;10;1'}");
    VParam v10("v10", p.get_child("a"));
    BOOST_CHECK (v10.size() == 10);
    BOOST_CHECK (v10.currentValue() == "1");
    v10.index++;
    BOOST_CHECK (v10.currentValue() == "2");
    v10.index++;
    BOOST_CHECK (v10.currentValue() == "3");

    p = ptreeFromJson("{'a': '1;2;0.2'}");
    VParam v11("v11", p.get_child("a"));
    BOOST_CHECK (v11.size() == 6);
    BOOST_CHECK (v11.currentValue() == "1");
    v11.index++;
    BOOST_CHECK (v11.currentValue() == "1.2");
    v11.index++;
    BOOST_CHECK (v11.currentValue() == "1.4");

    p = ptreeFromJson("{'a': '7;-2;-1.5'}");
    VParam v12("v12", p.get_child("a"));
    BOOST_CHECK (v12.size() == 7);
    BOOST_CHECK (v12.currentValue() == "7");
    v12.index++;
    BOOST_CHECK (v12.currentValue() == "5.5");
    v12.index++;
    BOOST_CHECK (v12.currentValue() == "4");
}

BOOST_AUTO_TEST_CASE ( test_vconfiguration )
{
    ptree p1 = ptreeFromJson("{'a': 1, 'b': 2}");
    VConfiguration vc1(p1);
    BOOST_CHECK (vc1.hasNext() == true);
    BOOST_CHECK (vc1.numberOfVariants() == 1);
    BOOST_CHECK (vc1.getNext() == p1);
    BOOST_CHECK (vc1.hasNext() == false);

    ptree p2 = ptreeFromJson("{'a': 1, 'b': 2, 'changing': 'a'}");
    VConfiguration vc2(p2);
    BOOST_CHECK_THROW (vc2.hasNext(), ConfigurationException);

    ptree p3 = ptreeFromJson("{'a': 1, 'b': 2, 'changing': ['a']}");
    VConfiguration vc3(p3);
    BOOST_CHECK_THROW (vc3.hasNext(), ConfigurationException);

    ptree p4 = ptreeFromJson("{'a': 1, 'b': 2, 'changing': ['blub']}");
    VConfiguration vc4(p4);
    BOOST_CHECK_THROW (vc4.hasNext(), ConfigurationException);

    ptree p5 = ptreeFromJson("{'a': [1], 'b': 2, 'changing': ['a']}");
    VConfiguration vc5(p5);
    BOOST_CHECK (vc5.hasNext() == true);
    BOOST_CHECK (vc5.numberOfVariants() == 1);
    BOOST_CHECK (vc5.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2, 'changing': ['a'], 'original': {'a': [1]}, 'changed': ['a']}"));
    BOOST_CHECK (vc5.hasNext() == false);

    ptree p6 = ptreeFromJson("{'a': [1,2,3], 'b': 2, 'changing': ['a']}");
    VConfiguration vc6(p6);
    BOOST_CHECK (vc6.hasNext() == true);
    BOOST_CHECK (vc6.numberOfVariants() == 3);
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2, 'changing': ['a'], 'original': {'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 2, 'changing': ['a'], 'original': {'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 2, 'changing': ['a'], 'original': {'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc6.hasNext() == false);

    ptree p7 = ptreeFromJson("{'a': [1,2,3], 'b': [1,2], 'changing': ['b','a']}");
    VConfiguration vc7(p7);
    BOOST_CHECK (vc7.hasNext() == true);
    BOOST_CHECK (vc7.numberOfVariants() == 6);
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 1, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['b','a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 1, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 1, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['b','a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 2, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 2, 'changing': ['b','a'], "
                "'original': {'b': [1,2], 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc7.hasNext() == false);

    ptree p8 = ptreeFromJson("{'a': [1,2,3], 'b': {'c': [1,2], 'd': 10}, 'changing': ['b.c','a']}");
    VConfiguration vc8(p8);
    BOOST_CHECK (vc8.hasNext() == true);
    BOOST_CHECK (vc8.numberOfVariants() == 6);
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 1, 'b': {'c': 1, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['b.c','a']}"));
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 2, 'b': {'c': 1, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 3, 'b': {'c': 1, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 1, 'b': {'c': 2, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['b.c','a']}"));
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 2, 'b': {'c': 2, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['a']}"));
    BOOST_CHECK (vc8.getNext() == ptreeFromJson(
                "{'a': 3, 'b': {'c': 2, 'd': 10}, 'changing': ['b.c','a'], "
                "'original': {'b': {'c': [1,2]}, 'a': [1,2,3]}, 'changed': ['a']}"));
}

BOOST_AUTO_TEST_CASE ( test_configurator )
{
    Configurator c1(jsonify("{'a': 1, 'b': 2}"));
    BOOST_CHECK (c1.numberOfConfigs() == 1);
    BOOST_CHECK (c1.hasNext() == true);
    BOOST_CHECK (c1.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2}"));
    BOOST_CHECK (c1.hasNext() == false);
    BOOST_CHECK_THROW (c1.getNext(), ConfigurationException);
    
    Configurator c2(jsonify("[{'a': 1, 'b': 2},{'c': 3},{'blah': 'blub'}]"));
    BOOST_CHECK (c2.numberOfConfigs() == 3);
    BOOST_CHECK (c2.hasNext() == true);
    BOOST_CHECK (c2.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2}"));
    BOOST_CHECK (c2.hasNext() == true);
    BOOST_CHECK (c2.getNext() == ptreeFromJson(
                "{'c': 3}"));
    BOOST_CHECK (c2.hasNext() == true);
    BOOST_CHECK (c2.getNext() == ptreeFromJson(
                "{'blah': 'blub'}"));
    BOOST_CHECK (c2.hasNext() == false);
    BOOST_CHECK_THROW (c2.getNext(), ConfigurationException);

    Configurator c3(jsonify("[{'a': [1,2,3], 'b': 1, 'changing': ['a']}, "
                            " {'a': 2, 'b': [1,2], 'c': [1,2], 'changing': ['b','c']}]"));
    BOOST_CHECK (c3.numberOfConfigs() == 7);
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 1, 'changing': ['a'], 'original': {'a': [1,2,3]}, "
                "'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 1, 'changing': ['a'], 'original': {'a': [1,2,3]}, "
                "'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 1, 'changing': ['a'], 'original': {'a': [1,2,3]}, "
                "'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 1, 'c': 1, 'changing': ['b','c'], "
                "'original': {'b': [1,2], 'c': [1,2]}, 'changed': ['b','c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 1, 'c': 2, 'changing': ['b','c'], "
                "'original': {'b': [1,2], 'c': [1,2]}, 'changed': ['c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 2, 'c': 1, 'changing': ['b','c'], "
                "'original': {'b': [1,2], 'c': [1,2]}, 'changed': ['b','c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 2, 'c': 2, 'changing': ['b','c'], "
                "'original': {'b': [1,2], 'c': [1,2]}, 'changed': ['c']}"));
    BOOST_CHECK (c3.hasNext() == false);
    BOOST_CHECK_THROW (c3.getNext(), ConfigurationException);
}

BOOST_AUTO_TEST_CASE ( test_store )
{
    Configurator c1(jsonify(
        "[{'a': 1, 'b': [2,3,4,5], 'c': 'blah', 'd': {'d1': 2.1, 'd2': 3.4}, 'changing': ['b']}, "
        " {'a': 2, 'b': [2,3,4,5,6,7], 'c': 'blah', 'd': {'d1': [2.1,3,4.1,6.4], 'd2': 3.4}, 'changing': ['b','d.d1']}, "
        " {'a': 3, 'b': '1;22;1', 'c': 'blah', 'd': {'d1': 2.1, 'd2': 3.4}, 'changing': ['b']}]"));
    rtree r1;
    r1.put("P.0", 0.8);
    r1.put("P.1", 0.7);
    r1.put("P.2", 0.6);
    r1.put("N.0.total", 2.0);
    r1.put("N.0.0", 0.5);
    r1.put("N.0.1", 0.5);
    r1.put("N.0.2", 0.5);
    r1.put("N.0.3", 0.5);
    r1.put("N.1.total", 2.0);
    r1.put("N.1.0", 0.4);
    r1.put("N.1.1", 0.6);
    r1.put("N.1.2", 0.4);
    r1.put("N.1.3", 0.6);
    r1.put("N.2.total", 2.0);
    r1.put("N.2.0", 0.3);
    r1.put("N.2.1", 0.7);
    r1.put("N.2.2", 0.3);
    r1.put("N.2.3", 0.7);
    Store o1;
    for (int i=0; i<50; i++)
        o1.store(c1.getNext(), r1);

    // real world example
    CQca s1;
    Configurator c2(jsonify(
        "{'model': 'bond', 't': 1, 'td': 0.2, 'V0': 100, " 
        "'beta': 1000, 'layout': {'type': 'wire', 'cells': 2, "
        "'a': '0.01;0.1;0.01', 'b': [0.02,0.03,0.04], 'Pext': 0.1, "
        "'epc': 2}, 'observables': {'P': 'all', 'N': 0}, "
        "'changing': ['layout.b', 'layout.a']}"));
    while (c2.hasNext())
    {
        const ptree& cc = c2.getNext();
        s1.setConfig(cc);
        rtree r = s1.measure();
        // TODO: we should be able to read the whole configuration back
        //ptree p = s1.getConfig();
        //o1.store(p,r);
        o1.store(cc,r);
    }
}
