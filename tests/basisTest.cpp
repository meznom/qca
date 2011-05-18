#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE basis test
#include <boost/test/unit_test.hpp>

#include "basis.hpp"
#include "system.hpp"

BOOST_AUTO_TEST_CASE ( construct_state_from_string )
{
    State s("10010");
    BOOST_CHECK (s.size() == 5);
    BOOST_CHECK (s.toString() == "10010");
    s = State("110");
    BOOST_CHECK (s.size() == 3);
    BOOST_CHECK (s.toString() == "110");

    BOOST_CHECK_THROW (State("100a1"), FermionicStateException);
}

BOOST_AUTO_TEST_CASE ( construct_simple_basis )
{
    Basis<Filter::SelectAll, Sorter::DontSort> b(4, Filter::SelectAll(), Sorter::DontSort());
    b.construct();
    BOOST_CHECK (b.size() == 16);

}

BOOST_AUTO_TEST_CASE ( costruct_basis_with_symmetry_operators )
{
    Basis<Filter::SelectAll, Sorter::DontSort> b1(4, Filter::SelectAll(), Sorter::DontSort());
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b1.addSymmetryOperator(N);
    b1.addSymmetryOperator(S);
    b1.construct();
    
    BOOST_CHECK (b1.size() == 16);
   
    char eArray1[][10] = {"0000", "0100", "0001", "1000", "0010", "0101", "1100", "0110", "1001", "0011", "1010", "1101", "0111", "1110", "1011", "1111"};
    std::vector<std::string> expected1(eArray1, eArray1+16);
    for (size_t i=0; i<b1.size(); i++)
        BOOST_CHECK (expected1[i] == b1(i).toString());

    Basis<Filter::SelectAll, Sorter::DontSort> b2(4, Filter::SelectAll(), Sorter::DontSort());
    b2.addSymmetryOperator(S);
    b2.addSymmetryOperator(N);
    b2.construct();

    BOOST_CHECK (b2.size() == 16);
   
    char eArray2[][10] = {"0101", "0100", "0001", "1101", "0111", "0000", "1100", "0110", "1001", "0011", "1111", "1000", "0010", "1110", "1011", "1010"};
    std::vector<std::string> expected2(eArray2, eArray2+16);
    for (size_t i=0; i<b2.size(); i++)
        BOOST_CHECK (expected2[i] == b2(i).toString());
    
    delete N;
    delete S;
}

BOOST_AUTO_TEST_CASE ( test_basis_blocks )
{
    Basis<Filter::SelectAll, Sorter::DontSort> b(8, Filter::SelectAll(), Sorter::DontSort());
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b.addSymmetryOperator(N);
    b.addSymmetryOperator(S);
    b.construct();

    b.sectorRange(4, 1);

    delete N;
    delete S;
}
