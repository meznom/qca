/*
 * g++ -lboost_unit_test_framework -o utilitiesTest utilitiesTest.cpp
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE utilities test
#include <boost/test/unit_test.hpp>

#include "../utilities.hpp"

BOOST_AUTO_TEST_CASE ( empty_description_item)
{
    DescriptionItem i;
    BOOST_CHECK (i.getValue() == "");

    BOOST_CHECK_THROW (static_cast<double>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<int>(i), ConversionException);
    BOOST_CHECK_THROW (static_cast<size_t>(i), ConversionException);
    BOOST_CHECK (i == "");
    BOOST_CHECK (static_cast<bool>(i) == false);
    BOOST_CHECK (i.isSet() == false);
    std::vector<int> v = i;
    BOOST_CHECK (v.size() == 0);
}
