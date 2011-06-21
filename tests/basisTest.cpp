#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE basis test
#include <boost/test/unit_test.hpp>

#define STORAGE_TYPE_OF_FERMIONIC_STATE uint32_t
#include "basis.hpp"

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
    Basis b(4);
    b.construct();
    BOOST_CHECK (b.size() == 16);
}

BOOST_AUTO_TEST_CASE ( costruct_basis_with_symmetry_operators )
{
    Basis b1(4);
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

    Basis b2(4);
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

BOOST_AUTO_TEST_CASE ( test_sectors_and_ranges_1 )
{
    Basis b(4);
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b.addSymmetryOperator(N);
    b.addSymmetryOperator(S);
    b.construct();

    BOOST_CHECK (b.getRanges().size() == 9);

    //std::vector<Range> rs = b.getRanges();
    //for (size_t i=0; i<rs.size(); i++)
    //    std::cerr << "NS> " << i << "  " << rs[i].b - rs[i].a << std::endl;

    Range r, er;
    r = b.getRangeOfSector(2, 0);
    er = Range(6, 10);
    BOOST_CHECK (r == er);

    for (size_t i=6; i<10; i++)
    {
        const State& s = b(i);
        BOOST_CHECK((*N)(s) == 2);
        BOOST_CHECK((*S)(s) == 0);
    }
    
    r = b.getRangeOfSector(3, 1);
    er = Range(13, 15);
    BOOST_CHECK (r == er);

    r = b.getRangeOfSector(2, 1); //empty/undefined sector
    er = Range(0, 0);
    BOOST_CHECK (r == er);

    delete N;
    delete S;
}

BOOST_AUTO_TEST_CASE ( test_sectors_and_ranges_2 )
{
    Basis b(4);
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b.addSymmetryOperator(S);
    b.addSymmetryOperator(N);
    b.construct();

    BOOST_CHECK (b.getRanges().size() == 9);

    //std::vector<Range> rs = b.getRanges();
    //for (size_t i=0; i<rs.size(); i++)
    //    std::cerr << "SN> " << i << "  " << rs[i].b - rs[i].a << std::endl;

    Range r, er;
    r = b.getRangeOfSector(0, 2);
    er = Range(6, 10);
    BOOST_CHECK (r == er);

    for (size_t i=6; i<10; i++)
    {
        const State& s = b(i);
        BOOST_CHECK((*N)(s) == 2);
        BOOST_CHECK((*S)(s) == 0);
    }
    
    r = b.getRangeOfSector(2, 2);
    er = Range(15, 16);
    BOOST_CHECK (r == er);

    r = b.getRangeOfSector(1, 2); //empty/undefined sector
    er = Range(0, 0);
    BOOST_CHECK (r == er);

    delete N;
    delete S;
}

BOOST_AUTO_TEST_CASE ( test_sectors_and_ranges_3 )
{
    Basis b(8);
    ParticleNumberSymmetryOperator* N = new ParticleNumberSymmetryOperator();
    SpinSymmetryOperator* S = new SpinSymmetryOperator();
    b.addSymmetryOperator(N);
    b.addSymmetryOperator(S);
    b.construct();

    Range r;
    r = b.getRangeOfSector(4, 4);
    BOOST_CHECK (r.b-r.a == 1);

    r = b.getRangeOfSector(4, 2);
    BOOST_CHECK (r.b-r.a == 16);

    r = b.getRangeOfSector(6, 2);
    BOOST_CHECK (r.b-r.a == 6);

    delete N;
    delete S;
}

BOOST_AUTO_TEST_CASE ( test_construct_basis_with_filter )
{
    //TODO
}

BOOST_AUTO_TEST_CASE ( performance_construct_basis )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    Basis b(22);
    ParticleNumberSymmetryOperator N;
    SpinSymmetryOperator S;
    b.addSymmetryOperator(&N);
    b.addSymmetryOperator(&S);
    
    startCPUTime = std::clock();
    b.construct();
    endCPUTime = std::clock();

    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction: " << cpuTime << "s" << std::endl;
}
