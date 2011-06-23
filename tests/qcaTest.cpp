#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qca test
#include <boost/test/unit_test.hpp>

#include <ctime>

#define STORAGE_TYPE_OF_FERMIONIC_STATE uint32_t
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

BOOST_AUTO_TEST_CASE ( test_construct_qca_quarterfilled_system )
{
    QcaQuarterFilling<ExternalPlain> s1(1);
    s1.construct();
    BOOST_CHECK (s1.basis.size() == 28);
    BOOST_CHECK (s1.basis.getRanges().size() == 3);

    QcaQuarterFilling<ExternalPlain> s2(2);
    s2.construct();
    BOOST_CHECK (s2.basis.size() == 28*28);
    BOOST_CHECK (s2.basis.getRanges().size() == 5);

#ifdef NDEBUG
    QcaQuarterFilling<ExternalPlain> s3(3);
    s3.construct();
    BOOST_CHECK (s3.basis.size() == 28*28*28);
    BOOST_CHECK (s3.basis.getRanges().size() == 7);
#endif
}

BOOST_AUTO_TEST_CASE ( test_qca_quarterfilling_system_for_some_parameters )
{
    //see tests/NotesAndTests.txt
    QcaQuarterFilling<ExternalDeadPlaquet> s1(1);
    s1.H.Vext = 0;
    s1.H.Pext = 0;
    s1.H.t = 1;
    s1.H.td = 0.2;
    s1.H.a = 1.0/250.0;
    s1.H.b = 1.75*s1.H.a;
    s1.H.V0 = 10.0 / s1.H.a;
    
    s1.H.Pext = 0;
    s1.construct();
    s1.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s1.ensembleAverage(1000, s1.P(0)), 0));

    s1.H.Pext = 0.1;
    s1.construct();
    s1.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s1.ensembleAverage(1000, s1.P(0)), 0.895, 0.001));

    QcaQuarterFilling<ExternalDeadPlaquet> s2(2);
    s2.H.Vext = 0;
    s2.H.Pext = 0;
    s2.H.t = 1;
    s2.H.td = 0.2;
    s2.H.a = 1.0/250.0;
    s2.H.b = 1.75*s2.H.a;
    s2.H.V0 = 10.0 / s2.H.a;
    
    s2.H.Pext = 0;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000, s2.P(0)), 0));
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000, s2.P(1)), 0));

    s2.H.Pext = 0.01;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(0)), 0.662, 0.001));
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(1)), 0.013, 0.001));

    s2.H.Pext = 0.2;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(0)), 0.998, 0.001));
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(1)), 0.020, 0.001));

    s2.H.Vext = 0;
    s2.H.Pext = 0;
    s2.H.t = 1;
    s2.H.td = 0;
    s2.H.a = 1.0/100.0;
    s2.H.b = 3*s2.H.a;
    s2.H.V0 = 10.0 / s2.H.a;
    
    s2.H.Pext = 0.1;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(0)), 0.414, 0.001));
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(1)), 0.379, 0.001));

    s2.H.Pext = -0.8;
    s2.construct();
    s2.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(0)), -0.944, 0.001));
    BOOST_CHECK (epsilonEqual(s2.ensembleAverage(1000000, s2.P(1)), -0.850, 0.001));
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_grandcanonical_system )
{
    QcaGrandCanonical<ExternalPlain> s1(1);
    s1.construct();
    BOOST_CHECK (s1.basis.size() == 256);
    BOOST_CHECK (s1.basis.getRanges().size() == 1+2+3+4+5+4+3+2+1);

    QcaGrandCanonical<ExternalPlain> s2(2);
    s2.construct();
    BOOST_CHECK (s2.basis.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( test_diagonalization_of_qca_grand_canonical_system )
{
    QcaGrandCanonical<ExternalPlain> s1(1);
    s1.construct();
    s1.H.diagonalize();
    DVector ev1_ = s1.H.eigenvalues;
    std::vector<double> ev1(ev1_.size());
    for (int i=0; i<ev1_.size(); i++)
        ev1[i] = ev1_(i);
    std::sort(ev1.begin(), ev1.end());
    
    s1.H.diagonalizeBlockWise();
    DVector ev2_ = s1.H.eigenvalues;
    std::vector<double> ev2(ev2_.size());
    for (int i=0; i<ev2_.size(); i++)
        ev2[i] = ev2_(i);
    std::sort(ev2.begin(), ev2.end());
    
    BOOST_CHECK (ev1.size() == ev2.size());
    for (size_t i=0; i<ev1.size(); i++)
        BOOST_CHECK (epsilonEqual(ev1[i], ev2[i]));

    QcaGrandCanonical<ExternalPlain> s2(2);
    s2.construct();

    /*
     * useful output
     */
    //const std::vector<Range>& rs = s2.basis.getRanges();
    //std::cerr << "For the two plaquet grand canonical QCA system we have " 
    //          << rs.size() << " ranges." << std::endl;
    //for (size_t i=0; i<rs.size(); i++)
    //    std::cerr << "   Size of Range " << i << ": " << rs[i].b-rs[i].a << std::endl;
    //std::cerr << std::endl;
    
    /*
     * doesn't work, because the dense matrix gets too big
     */
    //s2.H.diagonalize();

    /*
     * this takes a _long_ time, but works
     */
    //s2.H.diagonalizeBlockWise();
    //BOOST_CHECK (s2.H.eigenvalues.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( simple_sanity_checks_for_qca_grand_canonical_system )
{
    QcaGrandCanonical<ExternalDeadPlaquet> s1(1);
    s1.H.Vext = 0;
    s1.H.Pext = 0;
    s1.H.t = 1;
    s1.H.td = 0;
    s1.H.a = 1.0/160.0;
    s1.H.b = 3*s1.H.a;
    s1.H.V0 = 1000;
    s1.H.mu = -300;
    
    s1.H.Pext = 0;
    s1.construct();
    s1.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s1.ensembleAverage(1000000, s1.P(0)), 0));

    s1.H.Pext = 0.1;
    s1.construct();
    s1.H.diagonalizeBlockWise();
    BOOST_CHECK (epsilonEqual(s1.ensembleAverage(1000000, s1.N(0)), 2));
    double P = s1.ensembleAverage(100000, s1.P(0));
    BOOST_CHECK (P > 0.1 && P < 1);
}

BOOST_AUTO_TEST_CASE ( compare_performance_qca_quarterfilling_with_and_without_using_symmetries )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    QcaQuarterFilling<ExternalDeadPlaquet> s(2);

    startCPUTime = std::clock();
    s.construct();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalize();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalize of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeBlockWise();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeBlockWise of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;
}
