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
    //BOOST_CHECK (s1.measure() == rtree());

    CQca s2;
    s2.setConfig(ptreeFromJson("{}")); // default configuration
    BOOST_CHECK (s2.getConfig() == ptreeFromJson(
        "{'model': 'grandcanonical', 't': 1, 'td': 0, 'Vext': 0, 'V0': 1000, "
        "'mu': 0, 'epsilonr': 1, 'lambdaD': 0, 'q': 0, 'beta': 1, "
        "'layout': {'type': 'wire', 'cells': 1, 'a': 1, "
        "'b': 3, 'Pext': 0, 'epc': 2}, 'observables': {}}"));
    //std::cerr << ptreeToJson(s1.getConfig()) << std::endl;

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
    //std::cerr << ptreeToJson(s3.getConfig()) << std::endl;

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
                "{'a': 1, 'b': 2, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (vc5.hasNext() == false);

    ptree p6 = ptreeFromJson("{'a': [1,2,3], 'b': 2, 'changing': ['a']}");
    VConfiguration vc6(p6);
    BOOST_CHECK (vc6.hasNext() == true);
    BOOST_CHECK (vc6.numberOfVariants() == 3);
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 2, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (vc6.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 2, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (vc6.hasNext() == false);

    ptree p7 = ptreeFromJson("{'a': [1,2,3], 'b': [1,2], 'changing': ['b','a']}");
    VConfiguration vc7(p7);
    BOOST_CHECK (vc7.hasNext() == true);
    BOOST_CHECK (vc7.numberOfVariants() == 6);
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 1, 'changing': ['b','a'], 'changed': ['b','a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 1, 'changing': ['b','a'], 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 1, 'changing': ['b','a'], 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 1, 'b': 2, 'changing': ['b','a'], 'changed': ['b','a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 2, 'changing': ['b','a'], 'changed': ['a']}"));
    BOOST_CHECK (vc7.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 2, 'changing': ['b','a'], 'changed': ['a']}"));
    BOOST_CHECK (vc7.hasNext() == false);
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
                "{'a': 1, 'b': 1, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                "{'a': 2, 'b': 1, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                "{'a': 3, 'b': 1, 'changing': ['a'], 'changed': ['a']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 1, 'c': 1, 'changing': ['b','c'], 'changed': ['b','c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 1, 'c': 2, 'changing': ['b','c'], 'changed': ['c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 2, 'c': 1, 'changing': ['b','c'], 'changed': ['b','c']}"));
    BOOST_CHECK (c3.hasNext() == true);
    BOOST_CHECK (c3.getNext() == ptreeFromJson(
                " {'a': 2, 'b': 2, 'c': 2, 'changing': ['b','c'], 'changed': ['c']}"));
    BOOST_CHECK (c3.hasNext() == false);
    BOOST_CHECK_THROW (c3.getNext(), ConfigurationException);
}
