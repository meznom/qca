#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qca test
#include <boost/test/unit_test.hpp>

#include "qca.hpp"

bool epsilonEqual (double v, double w, double epsilon = 10E-10)
{
    return fabs(v-w) < epsilon;
}

BOOST_AUTO_TEST_CASE ( test_particle_number_per_plaquet_symmetry_operator )
{
    ParticleNumberPerPlaquetSymmetryOperator PP(4);
    State s(12);
   
    s = State("000000000000");
    BOOST_CHECK (PP(s) == 0);
    
    s = State("010000000000");
    BOOST_CHECK (PP(s) == 1);

    s = State("010000000011");
    BOOST_CHECK (PP(s) == 201);

    s = State("000010110010");
    BOOST_CHECK (PP(s) == 130);

    PP = ParticleNumberPerPlaquetSymmetryOperator();

    s = State("00000000000000000000000000000000");
    BOOST_CHECK (PP(s) == 0);
    
    s = State("00100000000000000000001101111111");
    BOOST_CHECK (PP(s) == 7201);

    s = State("00000000100000000000000000000000");
    BOOST_CHECK (PP(s) == 10);

    s = State("00000000000000110000000100000000");
    BOOST_CHECK (PP(s) == 120);
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_bond_system )
{
    QcaBond<ExternalPlain> s1(1);
    s1.construct();
    BOOST_CHECK (s1.basis.size() == 6);

    QcaBond<ExternalPlain> s2(2);
    s2.construct();
    BOOST_CHECK (s2.basis.size() == 36);

    QcaBond<ExternalPlain> s3(3);
    s3.construct();
    BOOST_CHECK (s3.basis.size() == 216);
}

BOOST_AUTO_TEST_CASE ( test_qca_bond_system_for_some_parameters )
{
    /*
     * These values here are coming from tests/NotesAndTests.txt.
     */
    QcaBond<ExternalPlain> s1(1);
    s1.H.Vext = 0;
    s1.H.t = 0.2;
    s1.H.td = 0.04;
    s1.H.a = 1;
    s1.H.b = 1;
    s1.H.V0 = 10;
    s1.construct();
    s1.H.diagonalize();
    BOOST_CHECK (s1.H.eigenvalues.size() == 6);
    double eArray1[6] = {0.43, 0.43, 0.92, 1.08 , 1.28 , 1.28};
    std::vector<double> expected1(eArray1, eArray1+6);
    for (size_t i=0; i<expected1.size(); i++)
        BOOST_CHECK (epsilonEqual(expected1[i], s1.H.eigenvalues(i), 0.01));

    QcaBond<ExternalPlain> s2(2);
    s2.H.Vext = 0.1;
    s2.H.t = 1;
    s2.H.td = 0;
    s2.H.a = 0.001;
    s2.H.b = 0.00175;
    s2.H.V0 = 1000000;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (s2.H.eigenvalues.size() == 36);
    double eArray2[5] = {2895.12, 2895.32, 2934.7, 2934.87, 2936.78};
    std::vector<double> expected2(eArray2, eArray2+5);
    for (size_t i=0; i<expected2.size(); i++)
        BOOST_CHECK (epsilonEqual(expected2[i], s2.H.eigenvalues(i), 0.01));

    QcaBond<ExternalDeadPlaquet> s3(1);
    s3.H.Vext = 0;
    s3.H.Pext = 0;
    s3.H.t = 1;
    s3.H.td = 0.2;
    s3.H.a = 1.0/250.0;
    s3.H.b = 1.75 * s3.H.a;
    s3.H.V0 = 10.0 / s3.H.a;

    s3.H.Pext = 0.01;
    s3.construct();
    s3.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s3.ensembleAverage(1000000, s3.P(0)), 0.29, 0.01));

    s3.H.Pext = 0.1;
    s3.construct();
    s3.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s3.ensembleAverage(1000000, s3.P(0)), 0.93, 0.01));
}
