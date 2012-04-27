#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cqca test
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/json_parser.hpp>
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

BOOST_AUTO_TEST_CASE ( test_configurable_qca_systems )
{
    std::string json1 = 
        "{                                      "
        "   \"t\": 1,                           "
        "   \"V0\": 10000,                      "
        "   \"layout\": {                       "
        "       \"type\": \"nonuniformwire\",   "
        "       \"a\": 10,                      "
        "       \"cells\": 2,                   "
        "       \"bs\": [1,2]                   "
        "   },                                  "
        "   \"observables\": {                  "
        "       \"P\": \"all\",                 "
        "       \"P2\": [0,1],                  "
        "       \"N\": 0,                       "
        "       \"E\": \"yes\"                  "
        "   }                                   "
        "}                                      ";
    ptree c;
    std::stringstream ss(json1);

    read_json(ss, c);
    for (ptree::const_iterator i = c.begin(); i!=c.end(); i++)
    {
        std::cerr << i->first << "  " << i->second.data() << std::endl;
    }

    CQca<QcaBond> s1;
    s1.setConfig(c);

    ptree c2 = s1.getConfig();
    std::stringstream ss2;
    write_json(ss2, c2);
    std::cerr << ss2.str() << std::endl;

    rtree r1 = s1.measure();
    //std::stringstream ss3;
    //write_json(ss3, r1);
    //std::cerr << ss3.str() << std::endl;
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
