#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE cqca test
#include <boost/test/unit_test.hpp>
#include <boost/property_tree/json_parser.hpp>
#include "cqca.hpp"

using boost::property_tree::ptree;
using namespace boost::property_tree::json_parser;

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

    ptree r1 = s1.measure();
    std::stringstream ss3;
    write_json(ss3, r1);
    std::cerr << ss3.str() << std::endl;
}
