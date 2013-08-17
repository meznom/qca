#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qca test
#include <boost/test/unit_test.hpp>

#define STORAGE_TYPE_OF_FERMIONIC_STATE uint32_t
#include "qca.hpp"

BOOST_AUTO_TEST_CASE ( test_construct_qca_ising_system )
{
    QcaIsing s;
    s.l.wire(1,1,1,1);
    s.update();
}
