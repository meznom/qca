#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE basis test
#include <boost/test/unit_test.hpp>

#include "basis.hpp"
#include "system.hpp"

BOOST_AUTO_TEST_CASE ( construct_simple_basis )
{
    Basis<Filter::SelectAll, Sorter::DontSort> b(4, Filter::SelectAll(), Sorter::DontSort());
    b.construct();
    BOOST_CHECK (b.size() == 16);

}

BOOST_AUTO_TEST_CASE ( costruct_with_symmetry_operators )
{
    Basis<Filter::SelectAll, Sorter::DontSort> b(4, Filter::SelectAll(), Sorter::DontSort());
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b.addSymmetryOperator(N);
    b.addSymmetryOperator(S);
    b.construct();
    BOOST_CHECK (b.size() == 16);
    b.dump();
    delete N;
    delete S;
}
