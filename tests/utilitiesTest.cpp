/*
 * g++ -lboost_unit_test_framework -o utilitiesTest utilitiesTest.cpp
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE utilities test
#include <boost/test/unit_test.hpp>

#include "utilities.hpp"

BOOST_AUTO_TEST_CASE ( empty_DescriptionItem )
{
    DescriptionItem i;
    BOOST_CHECK (i.getValue() == "");

    BOOST_CHECK_THROW (static_cast<double>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<int>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);
    BOOST_CHECK (i == "");
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);
    BOOST_CHECK (i.isSet() == false);
    std::vector<int> v = i;
    BOOST_CHECK (v.size() == 0);
}

BOOST_AUTO_TEST_CASE ( DescriptionItem_type_conversion )
{
    DescriptionItem i("4.2");

    BOOST_CHECK (i == "4.2");
    BOOST_CHECK (static_cast<double>(i) == 4.2);
    BOOST_CHECK_THROW (static_cast<int>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);
    std::vector<double> vd = i;
    BOOST_CHECK (vd.size() == 1);
    BOOST_CHECK (vd[0] == 4.2);
    BOOST_CHECK_THROW (std::vector<int> vi = i, ConversionException);

    i = -7.2;
    BOOST_CHECK (i == "-7.2");
    BOOST_CHECK (static_cast<double>(i) == -7.2);

    i = -5;
    BOOST_CHECK (i == "-5");
    BOOST_CHECK (static_cast<int>(i) == -5);
    BOOST_CHECK (static_cast<double>(i) == -5.0);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);

    i = 230;
    BOOST_CHECK (i == "230");
    BOOST_CHECK (static_cast<int>(i) == 230);
    BOOST_CHECK (static_cast<double>(i) == 230.0);
    BOOST_CHECK (static_cast<size_t>(i) == 230);

    i = "0";
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i = "";
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i = 4;
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);
    
    i = "false";
    BOOST_CHECK (static_cast<bool>(i) == false);
    BOOST_CHECK (static_cast<bool>(i) != true);

    i = "true";
    BOOST_CHECK (static_cast<bool>(i) == true);
    BOOST_CHECK (static_cast<bool>(i) != false);

    //TODO: test / check what happens when the number is too large for int or
    //size_t, etc
}

BOOST_AUTO_TEST_CASE ( DescriptionItem_default_value_and_isSet )
{
    DescriptionItem i;

    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK_THROW (static_cast<bool>(i), ConversionException);

    i.setDefault(7);
    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK (static_cast<int>(i) == 7);

    i.setDefault(true);
    BOOST_CHECK (i.isSet() == false);
    BOOST_CHECK (static_cast<bool>(i) == true);

    i.setDefault(-8.4);
    BOOST_CHECK (i.get(5) == 5);
    BOOST_CHECK (i.get<double>() == -8.4);
    BOOST_CHECK (i.isSet() == false);

    i = -8.4;
    BOOST_CHECK (i.get<double>(5) == -8.4);
    BOOST_CHECK (i.get<double>() == -8.4);
    BOOST_CHECK (i.isSet() == true);
}

BOOST_AUTO_TEST_CASE ( DescriptionItem_lists )
{
    DescriptionItem i;
    std::vector<double> v;

    i = "[]";
    v = i;
    BOOST_CHECK (v.size() == 0);

    i = "";
    v = i;
    BOOST_CHECK (v.size() == 0);

    i = "()";
    BOOST_CHECK_THROW (v = i, ConversionException);

    i = "[1]";
    v = i;
    BOOST_CHECK (v.size() == 1);
    BOOST_CHECK (v[0] == 1);

    i = "1";
    v = i;
    BOOST_CHECK (v.size() == 1);
    BOOST_CHECK (v[0] == 1);

    i = "[10,20,30,40,50]";
    v = i;
    BOOST_CHECK (v.size() == 5);
    BOOST_CHECK (v[3] == 40);

    /*
     * TODO: not yet implemented
     */
    /*
    i = "10,20,30,40,50";
    v = i;
    BOOST_CHECK (v.size() == 5);
    BOOST_CHECK (v[3] == 40);

    i = "[10,20,30,40,50,]";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,40,50,";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,40,,50";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10,20,30,4 0,50";
    BOOST_CHECK_THROW (v = i, ConversionException);

    i = "10, 20,30,40 ,50";
    v = i;
    BOOST_CHECK (v.size() == 5);
    BOOST_CHECK (v[3] == 40);

    i = "(10)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(10,20)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(10,20,30,40)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "10:20:30:40";
    BOOST_CHECK_THROW (v = i, ConversionException);

    i = "(2,1,0.1)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "(1,2,-0.1)";
    BOOST_CHECK_THROW (v = i, ConversionException);
    i = "1:2:-0.1";
    BOOST_CHECK_THROW (v = i, ConversionException);
    
    i = "(1,2.01,0.1)";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (v[2] == 1.2);
    
    i = "1:2.01:0.1";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (v[2] == 1.2);
    
    i = "(2,0.999,0.1)";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (v[2] == 1.8);
    
    i = "2:0.999:0.1";
    v = i;
    BOOST_CHECK (v.size() == 11);
    BOOST_CHECK (v[2] == 1.8);
    */
}
