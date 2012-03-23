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
    QcaBond s1(1, false);
    BOOST_CHECK (s1.basis.size() == 6);

    QcaBond s2(2, false);
    BOOST_CHECK (s2.basis.size() == 36);

    QcaBond s3(3, false);
    BOOST_CHECK (s3.basis.size() == 216);
}

BOOST_AUTO_TEST_CASE ( test_qca_bond_system_for_some_parameters )
{
    /*
     * These values here are coming from tests/NotesAndTests.txt.
     */
    QcaBond s1(1, false);
    s1.Vext = 0;
    s1.t = 0.2;
    s1.td = 0.04;
    s1.a = 1;
    s1.b = 1;
    s1.V0 = 10;
    s1.layout = Layout();
    s1.layout.addWireNoDriver(0,0, 1, s1.a, s1.b, 0);
    s1.H.construct();
    s1.H.diagonalizeNoSymmetries();
    BOOST_CHECK (s1.energies().size() == 6);
    double eArray1[6] = {0.43, 0.43, 0.92, 1.08 , 1.28 , 1.28};
    std::vector<double> expected1(eArray1, eArray1+6);
    for (size_t i=0; i<expected1.size(); i++)
        BOOST_CHECK (epsilonEqual(expected1[i], s1.energies()(i), 0.01));

    QcaBond s2(2, false);
    s2.Vext = 0.1;
    s2.t = 1;
    s2.td = 0;
    s2.a = 0.001;
    s2.b = 0.00175;
    s2.V0 = 1000000;
    s2.update();
    BOOST_CHECK (s2.energies().size() == 36);
    double eArray2[5] = {2895.12, 2895.32, 2934.7, 2934.87, 2936.78};
    std::vector<double> expected2(eArray2, eArray2+5);
    for (size_t i=0; i<expected2.size(); i++)
        BOOST_CHECK (epsilonEqual(expected2[i], s2.energies()(i), 0.01));

    QcaBond s3(1, true);
    s3.Vext = 0;
    s3.Pext = 0;
    s3.t = 1;
    s3.td = 0.2;
    s3.a = 1.0/250.0;
    s3.b = 1.75 * s3.a;
    s3.V0 = 10.0 / s3.a;

    s3.Pext = 0.01;
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000000, s3.P(0)), 0.29, 0.01));

    s3.Pext = 0.1;
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000000, s3.P(0)), 0.93, 0.01));

    QcaBond s4(3, true);
    s4.Vext = 0;
    s4.Pext = 1;
    s4.t = 1;
    s4.td = 0;
    s4.ti = 0;
    s4.a = 0.01;
    s4.b =  0.023;
    s4.q = 0;
    s4.epsilonr = QCA_NATURAL_EPSILON_R; // natural units
    s4.lambdaD = 0; 
    s4.V0 = 1000;
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(0)), 0.670243, 1E-5));
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(1)), 0.462795, 1E-5));
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(2)), 0.295182, 1E-5));
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_fixed_charge_system )
{
    QcaFixedCharge s1(1, 2, false);
    BOOST_CHECK (s1.basis.size() == 28);
    BOOST_CHECK (s1.basis.getRanges().size() == 3);

    QcaFixedCharge s2(2, 2, false);
    BOOST_CHECK (s2.basis.size() == 28*28);
    BOOST_CHECK (s2.basis.getRanges().size() == 5);

#ifdef NDEBUG
    QcaFixedCharge s3(3, 2, false);
    BOOST_CHECK (s3.basis.size() == 28*28*28);
    BOOST_CHECK (s3.basis.getRanges().size() == 7);
#endif

    QcaFixedCharge s4(1, 6, false);
    BOOST_CHECK (s4.basis.size() == 28);

    QcaFixedCharge s5(2, 6, false);
    BOOST_CHECK (s5.basis.size() == 28*28);
}

BOOST_AUTO_TEST_CASE ( test_qca_fixed_charge_system_for_some_parameters )
{
    //see tests/NotesAndTests.txt
    QcaFixedCharge s1(1, 2, true);
    s1.Vext = 0;
    s1.Pext = 0;
    s1.t = 1;
    s1.td = 0.2;
    s1.a = 1.0/250.0;
    s1.b = 1.75*s1.a;
    s1.V0 = 10.0 / s1.a;
    
    s1.Pext = 0;
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.N(0)), 2));
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.P(0)), 0));

    s1.Pext = 0.1;
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.P(0)), 0.895, 0.001));


    QcaFixedCharge s2(2, 2, true);
    s2.Vext = 0;
    s2.Pext = 0;
    s2.t = 1;
    s2.td = 0.2;
    s2.a = 1.0/250.0;
    s2.b = 1.75*s2.a;
    s2.V0 = 10.0 / s2.a;
    
    s2.Pext = 0;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.N(0)), 2));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.N(1)), 2));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.P(0)), 0));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.P(1)), 0));

    s2.Pext = 0.01;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.662, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.013, 0.001));

    s2.Pext = 0.2;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.998, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.020, 0.001));

    s2.Vext = 0;
    s2.Pext = 0;
    s2.t = 1;
    s2.td = 0;
    s2.a = 1.0/100.0;
    s2.b = 3*s2.a;
    s2.V0 = 10.0 / s2.a;
    
    s2.Pext = 0.1;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.414, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.379, 0.001));

    s2.Pext = -0.8;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), -0.944, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), -0.850, 0.001));


    QcaFixedCharge s3(1, 6, true);
    s3.Vext = 0;
    s3.Pext = 0;
    s3.t = 1;
    s3.td = 0.2;
    s3.a = 1.0/250.0;
    s3.b = 4*s3.a;
    s3.V0 = 10.0 / s3.a;
    
    s3.Pext = 0;
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.N(0)), 6));
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.P(0)), 0));

    s3.Pext = 0.1;
    s3.update();
    BOOST_CHECK (s3.measure(1000, s3.P(0)) > 0.1);
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.P(0)), 0.213, 0.001));


    QcaFixedCharge s4(2, 6, true);
    s4.Vext = 0;
    s4.Pext = 0;
    s4.t = 1;
    s4.td = 0.2;
    s4.a = 1.0/250.0;
    s4.b = 5*s4.a;
    s4.V0 = 10.0 / s4.a;

    s4.Pext = 0;
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.N(0)), 6));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(0)), 0, 10E-8));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(1)), 0, 10E-8));

    s4.Pext = 0.1;
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(0)), 0.268, 0.001));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(1)), 0.212, 0.001));
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_grandcanonical_system )
{
    QcaGrandCanonical s1(1, 2, false);
    BOOST_CHECK (s1.basis.size() == 256);
    BOOST_CHECK (s1.basis.getRanges().size() == 1+2+3+4+5+4+3+2+1);

    QcaGrandCanonical s2(2, 2, false);
    BOOST_CHECK (s2.basis.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( test_diagonalization_of_qca_grand_canonical_system )
{
    QcaGrandCanonical s1(1, 2, false);
    s1.H.construct();
    s1.H.diagonalizeNoSymmetries();
    DVector ev1_ = s1.energies();
    std::vector<double> ev1(ev1_.size());
    for (int i=0; i<ev1_.size(); i++)
        ev1[i] = ev1_(i);
    std::sort(ev1.begin(), ev1.end());
    
    s1.H.diagonalizeUsingSymmetries();
    DVector ev2_ = s1.energies();
    std::vector<double> ev2(ev2_.size());
    for (int i=0; i<ev2_.size(); i++)
        ev2[i] = ev2_(i);
    std::sort(ev2.begin(), ev2.end());

    s1.H.diagonalizeUsingSymmetriesBySectors();
    DVector ev3_ = s1.energies();
    std::vector<double> ev3(ev3_.size());
    for (int i=0; i<ev3_.size(); i++)
        ev3[i] = ev3_(i);
    std::sort(ev3.begin(), ev3.end());
    
    BOOST_CHECK (ev1.size() == ev2.size());
    BOOST_CHECK (ev1.size() == ev3.size());
    for (size_t i=0; i<ev1.size(); i++)
    {
        BOOST_CHECK (epsilonEqual(ev1[i], ev2[i]));
        BOOST_CHECK (epsilonEqual(ev1[i], ev3[i]));
    }

    QcaGrandCanonical s2(2, 2, false);
    s2.H.construct();

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
    //s2.H.diagonalizeNoSymmetries();

    /*
     * this takes a _long_ time, but works
     */
    //s2.H.diagonalizeUsingSymmetries();
    //BOOST_CHECK (s2.H.eigenvalues.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( simple_sanity_checks_for_qca_grand_canonical_system )
{
    QcaGrandCanonical s1(1, 2, true);
    s1.Vext = 0;
    s1.Pext = 0;
    s1.t = 1;
    s1.td = 0;
    s1.a = 1.0/160.0;
    s1.b = 3*s1.a;
    s1.V0 = 1000;
    s1.mu = 300;
    
    s1.Pext = 0;
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000000, s1.P(0)), 0));

    s1.Pext = 0.1;
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000000, s1.N(0)), 2));
    double P = s1.measure(100000, s1.P(0));
    BOOST_CHECK (P > 0.1 && P < 1);


    QcaGrandCanonical s2(1, 6, true);
    s2.Vext = 0;
    s2.Pext = 0;
    s2.t = 1;
    s2.td = 0;
    s2.a = 1.0/160.0;
    s2.b = 4*s2.a;
    s2.V0 = 1000;
    s2.mu = 1800;
    
    s2.Pext = 0;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0));

    s2.Pext = 1;
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.N(0)), 6));
    P = s2.measure(100000, s2.P(0));
    BOOST_CHECK (P > 0.1 && P < 1);
}

BOOST_AUTO_TEST_CASE ( compare_performance_qca_fixed_charge_with_and_without_using_symmetries )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;

    startCPUTime = std::clock();
    QcaFixedCharge s(2, 2, true);
    s.H.construct();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeNoSymmetries();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeNoSymmetries of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetries();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetries of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetriesBySectors();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetriesBySectors of two plaquet QCA quarterfilled system: " 
              << cpuTime << "s" << std::endl;
}

BOOST_AUTO_TEST_CASE ( test_scaling_of_parameters_for_grand_canonical_system )
{
    QcaGrandCanonical s1(1, 2, false);
    s1.Vext = 0;
    s1.Pext = 0;
    s1.t = 1;
    s1.td = 0;
    s1.a = 1.0/160.0;
    s1.b = 3*s1.a;
    s1.V0 = 1000;
    s1.mu = 300;
    s1.q = 0;

    s1.update();
    double N1 = s1.measure(0.1, s1.N(0));
    double P1 = s1.measure(0.1, s1.P(0));

    QcaGrandCanonical s2(1, 2, false);
    s2.Vext = 0;
    s2.Pext = 0;
    s2.t = 1 * 10;
    s2.td = 0;
    s2.a = 1.0/160.0 * 1.0/10.0;
    s2.b = 3*s2.a;
    s2.V0 = 1000 * 10;
    s2.mu = 300 * 10;
    s2.q = 0;

    s2.update();
    double N2 = s2.measure(0.1 * 1.0/10.0, s2.N(0));
    double P2 = s2.measure(0.1 * 1.0/10.0, s2.P(0));

    BOOST_CHECK (epsilonEqual(N1, N2));
    BOOST_CHECK (epsilonEqual(P1, P2));


    QcaGrandCanonical s3(1, 2, false);
    s3.Vext = 1;
    s3.Pext = 0;
    s3.t = 1;
    s3.td = 0;
    s3.a = 1.0/160.0;
    s3.b = 3*s3.a;
    s3.V0 = 1000;
    s3.mu = 300;
    s3.q = 0;

    s3.update();
    double N3 = s3.measure(0.1, s3.N(0));
    double P3 = s3.measure(0.1, s3.P(0));

    QcaGrandCanonical s4(1, 2, false);
    s4.Vext = 1 * 1000;
    s4.Pext = 0;
    s4.t = 1 * 1000;
    s4.td = 0;
    s4.a = 1.0/160.0 * 1.0/1000.0;
    s4.b = 3*s4.a;
    s4.V0 = 1000 * 1000;
    s4.mu = 300 * 1000;
    s4.q = 0;

    s4.update();
    double N4 = s4.measure(0.1 * 1.0/1000.0, s4.N(0));
    double P4 = s4.measure(0.1 * 1.0/1000.0, s4.P(0));

    BOOST_CHECK (epsilonEqual(N3, N4));
    BOOST_CHECK (epsilonEqual(P3, P4));
}

BOOST_AUTO_TEST_CASE ( test_compensation_charge )
{
    /*
     * Grand canonical, one plaquet: Should be four electrons per plaquet at 
     * mu = 0.
     */
    QcaGrandCanonical s1(1, 2, false);
    s1.Vext = 0;
    s1.Pext = 0;
    s1.t = 1;
    s1.td = 0;
    s1.a = 1.0/160.0;
    s1.b = 3*s1.a;
    s1.V0 = 0;
    s1.mu = 0;
    s1.q = 1;

    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(10000, s1.N(0)), 4));

    /*
     * Fixed, two plaquets: For the case that doubly occupied states are
     * completely gapped out, the 2e per plaquet system without compensation
     * charges and the 6e per plaquet system with compensation charges should
     * behave identical.
     */
    QcaFixedCharge s2(2, 2, false);
    s2.Vext = 0.1;
    s2.Pext = 0;
    s2.t = 1;
    s2.td = 0;
    s2.a = 1.0/100.0;
    s2.b = 2*s2.a;
    s2.V0 = 10000;
    s2.mu = 0;
    s2.q = 0;

    s2.update();
    double P21 = s2.measure(10, s2.P(0));
    double P22 = s2.measure(10, s2.P(1));

    QcaFixedCharge s3(2, 6, false);
    s3.Vext = 0.1;
    s3.Pext = 0;
    s3.t = 1;
    s3.td = 0;
    s3.a = 1.0/100.0;
    s3.b = 2*s3.a;
    s3.V0 = 10000;
    s3.mu = 0;
    s3.q = 1;

    s3.update();
    double P31 = s3.measure(10, s3.P(0));
    double P32 = s3.measure(10, s3.P(1));

    //std::cerr << "P21 = " << P21 << ", P22 = " << P22 << std::endl;
    //std::cerr << "P31 = " << P31 << ", P32 = " << P32 << std::endl;

    BOOST_CHECK (epsilonEqual(P21, P31));
    BOOST_CHECK (epsilonEqual(P22, P32));

    /*
     * Grand canonical, one plaquet: The two electron dead plaquet system
     * without compensation charge should be the same as the six electron dead
     * plaquet system with compensation charge, provided we set the chemical
     * potential correctly.
     */
    QcaGrandCanonical s4(1, 6, true);
    s4.Vext = 0;
    s4.Pext = 1;
    s4.t = 1;
    s4.td = 0;
    s4.a = 1.0/100.0;
    s4.b = 4*s4.a;
    s4.V0 = 1000;
    s4.mu = 1200;
    s4.q = 1;

    s4.update();
    double P41 = s4.measure(10, s4.P(0));

    QcaGrandCanonical s5(1, 2, true);
    s5.Vext = 0;
    s5.Pext = 1;
    s5.t = 1;
    s5.td = 0;
    s5.a = 1.0/100.0;
    s5.b = 4*s5.a;
    s5.V0 = 1000;
    s5.mu = 200; //1000 less than s4.mu (because V0=1000)
    s5.q = 0;

    s5.update();
    double P51 = s5.measure(10, s5.P(0));

    //std::cerr << "P41 = " << P41 << ", P51 = " << P51 << std::endl;

    BOOST_CHECK (epsilonEqual(P41, P51));

    /*
     * Fixed, two plaquets: The 2e and 6e dead plaquet systems should be
     * equivalent.
     */
    QcaFixedCharge s6(2, 2, true);
    s6.Vext = 0;
    s6.Pext = 1;
    s6.t = 1;
    s6.td = 0;
    s6.a = 1.0/100.0;
    s6.b = 4*s6.a;
    s6.V0 = 1000;
    s6.mu = 0;
    s6.q = 0;

    s6.update();
    double P61 = s6.measure(10, s6.P(0));
    double P62 = s6.measure(10, s6.P(1));

    QcaFixedCharge s7(2, 6, true);
    s7.Vext = 0;
    s7.Pext = 1;
    s7.t = 1;
    s7.td = 0;
    s7.a = 1.0/100.0;
    s7.b = 4*s7.a;
    s7.V0 = 1000;
    s7.mu = 0;
    s7.q = 1;

    s7.update();
    double P71 = s7.measure(10, s7.P(0));
    double P72 = s7.measure(10, s7.P(1));

    //std::cerr << "P61 = " << P61 << ", P62 = " << P62 << std::endl;
    //std::cerr << "P71 = " << P71 << ", P72 = " << P72 << std::endl;

    BOOST_CHECK (epsilonEqual(P61, P71));
    BOOST_CHECK (epsilonEqual(P62, P72));
}

BOOST_AUTO_TEST_CASE ( compare_polarization_and_polarization2 )
{
    QcaGrandCanonical s1(1, 6, true);
    s1.Vext = 0;
    s1.Pext = 1;
    s1.t = 1;
    s1.td = 0;
    s1.a = 1.0/100.0;
    s1.b = 4*s1.a;
    s1.V0 = 1000;
    s1.mu = 1200;
    s1.q = 1;

    s1.update();
    double P11 = s1.measure(10, s1.P(0));
    double P12 = s1.measurePolarization2(10, 0);

    //std::cerr << "P11 = " << P11 << ", P12 = " << P12 << std::endl;
    BOOST_CHECK (epsilonEqual(s1.measure(10, s1.N(0)), 6));
    BOOST_CHECK (epsilonEqual(P11, P12, 1e-5));

    QcaGrandCanonical s2(1, 6, true);
    s2.Vext = 0;
    s2.Pext = 1;
    s2.t = 1;
    s2.td = 0;
    s2.a = 1.0/100.0;
    s2.b = 4*s2.a;
    s2.V0 = 10;
    s2.mu = 200;
    s2.q = 1;

    s2.update();
    double P21 = s2.measure(10, s2.P(0));
    double P22 = s2.measurePolarization2(10, 0);

    //std::cerr << "P21 = " << P21 << ", P22 = " << P22 << std::endl;
    BOOST_CHECK (epsilonEqual(s2.measure(10, s2.N(0)), 6));
    BOOST_CHECK (P21 > 0.9);
    BOOST_CHECK (P22 < 0.1);

    QcaFixedCharge s3(2, 6, true);
    s3.Vext = 0;
    s3.Pext = 1;
    s3.t = 1;
    s3.td = 0;
    s3.a = 1.0/100.0;
    s3.b = 4*s3.a;
    s3.V0 = 1000;
    s3.mu = 0;
    s3.q = 1;

    s3.update();
    double P311 = s3.measure(10, s3.P(0));
    double P312 = s3.measure(10, s3.P(1));
    double P321= s3.measurePolarization2(10, 0);
    double P322= s3.measurePolarization2(10, 1);

    BOOST_CHECK (epsilonEqual(P311, P321, 1e-5));
    BOOST_CHECK (epsilonEqual(P312, P322, 1e-5));

    QcaFixedCharge s4(2, 6, true);
    s4.Vext = 0;
    s4.Pext = 1;
    s4.t = 1;
    s4.td = 0;
    s4.a = 1.0/100.0;
    s4.b = 4*s4.a;
    s4.V0 = 10;
    s4.mu = 0;
    s4.q = 1;

    s4.update();
    double P411 = s4.measure(10, s4.P(0));
    double P412 = s4.measure(10, s4.P(1));
    double P421= s4.measurePolarization2(10, 0);
    double P422= s4.measurePolarization2(10, 1);

    BOOST_CHECK (std::fabs(P411) > 0.9);
    BOOST_CHECK (std::fabs(P412) > 0.9);
    BOOST_CHECK (std::fabs(P421) < 0.1);
    BOOST_CHECK (std::fabs(P422) < 0.1);
}

BOOST_AUTO_TEST_CASE ( test_basic_layouts )
{
    Layout l1;
    l1.addDot(0,0)
      .addDot(0,1)
      .addDot(1,1)
      .addDot(1,0)
      .addCharge(-1,1,1);
    BOOST_CHECK (l1.r(0,1) == 1);
    BOOST_CHECK (l1.r(1,3) == std::sqrt(2));
    BOOST_CHECK (l1.r_charge_dot(0,3) == std::sqrt(5));

    Layout l2;
    l2.addCell(0,0,0.2)
      .addCell(0,2,0.2)
      .addDriverCell(1,0,0.2,1)
      .addDot(1,1)
      .addCharge(0,-1,0.5);
    BOOST_CHECK(epsilonEqual(l2.r(0,3), 0.2));
    BOOST_CHECK(epsilonEqual(l2.r(4,6), std::sqrt(0.2*0.2+0.2*0.2)));
    BOOST_CHECK(epsilonEqual(l2.r(0,4), 2));
    BOOST_CHECK(epsilonEqual(l2.r(1,7), std::sqrt(1.8*1.8+0.2*0.2)));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(0,0), 1));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(2,1), 1.2));
    BOOST_CHECK(epsilonEqual(l2.r(0,8), std::sqrt(2)));
    BOOST_CHECK(epsilonEqual(l2.r(6,8), std::sqrt(1.2*1.2+0.8*0.8)));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(4,3), std::sqrt(1+0.2*0.2)));
    BOOST_CHECK(l2.N_dots() == 9);
    BOOST_CHECK(l2.N_charges() == 5);

    Layout l3;
    l3.addWire(10,10, 4, 0.2, 0.4, 0.7);
    BOOST_CHECK(l3.N_dots() == 16);
    BOOST_CHECK(l3.N_charges() == 4);
    BOOST_CHECK(epsilonEqual(
                0.5*(l3.charge(0)+l3.charge(2)-l3.charge(1)-l3.charge(3)), 
                0.7
    ));
    BOOST_CHECK(epsilonEqual(
                l3.r(4,14), 
                std::sqrt((3*0.2+2*0.4)*(3*0.2+2*0.4)+0.2*0.2)
    ));
    BOOST_CHECK(epsilonEqual(
                l3.r_charge_dot(0,0), 
                0.2+0.4
    ));
    BOOST_CHECK(epsilonEqual(
                l3.r_charge_dot(2,4), 
                std::sqrt((2*0.4+1*0.2)*(2*0.4+1*0.2)+0.2*0.2)
    ));

    Layout l4;
    std::vector<double> bs;
    bs.push_back(0.4);
    bs.push_back(0.4);
    bs.push_back(0.4);
    bs.push_back(0.4);
    l4.addNonuniformWire(3.2,3.2, 4, 0.2, bs, -0.2, 6);
    BOOST_CHECK(l4.N_dots() == 16);
    BOOST_CHECK(l4.N_charges() == 4);
    BOOST_CHECK(epsilonEqual(
                0.5*(l4.charge(0)+l4.charge(2)-l4.charge(1)-l4.charge(3)), 
                -0.2
    ));
    BOOST_CHECK(epsilonEqual(
                l4.r(4,14), 
                std::sqrt((3*0.2+2*0.4)*(3*0.2+2*0.4)+0.2*0.2)
    ));
    BOOST_CHECK(epsilonEqual(
                l4.r_charge_dot(0,0), 
                0.2+0.4
    ));
    BOOST_CHECK(epsilonEqual(
                l4.r_charge_dot(2,4), 
                std::sqrt((2*0.4+1*0.2)*(2*0.4+1*0.2)+0.2*0.2)
    ));

    Layout l5;
    bs.clear();
    bs.push_back(1);
    bs.push_back(1);
    bs.push_back(1.2);
    l5.addNonuniformWire(0,0, 3, 0.1, bs, 0);
    BOOST_CHECK(epsilonEqual(
                l5.r_charge_dot(0,0),
                1.1
    ));
    BOOST_CHECK(epsilonEqual(
                l5.r(7,8),
                1.2
    ));
    BOOST_CHECK(epsilonEqual(
                l5.r(0,9),
                std::sqrt((2*0.1+1+1.2)*(2*0.1+1+1.2)+0.1*0.1)
    ));

    Layout l6;
    l6.addDriverCell(0,0, 1, 0.3, 6);
    BOOST_CHECK(epsilonEqual(
                0.5*(l6.charge(0)+l6.charge(2)-l6.charge(1)-l6.charge(3)), 
                0.3
    ));

    Layout l7;
    l7.addDriverCell(0,0, 1, -0.12, 6);
    BOOST_CHECK(epsilonEqual(
                0.5*(l7.charge(0)+l7.charge(2)-l7.charge(1)-l7.charge(3)), 
                -0.12
    ));
}
