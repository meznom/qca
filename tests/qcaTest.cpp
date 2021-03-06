/*
 * Copyright (c) 2011-2014 Burkhard Ritter
 * This code is distributed under the two-clause BSD License.
 */
#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qca test
#include <boost/test/unit_test.hpp>
#include <ctime>
#include "qca.hpp"

// TODO: replace Vext tests (commented for now) with tests that don't need Vext

class Wire2e : public Layout
{
public:
    Wire2e (int N_p, double a = 1, double b = 3, double P = 0)
    {
        epc = epc2;
        wire(N_p, a, b, P, epc2);
    }
};

class Wire6e : public Layout
{
public:
    Wire6e (int N_p, double a = 1, double b = 3, double P = 0)
    {
        epc = epc6;
        wire(N_p, a, b, P, epc6);
    }
};

class WireNoDriver2e : public Layout
{
public:
    WireNoDriver2e (int N_p, double a = 1, double b = 3)
    {
        epc = epc2;

        double r_x = 0;
        double r_y = 0;
        
        for (int i=0; i<N_p; i++)
            addCell(r_x+i*(a+b), r_y, a);
    }
};

class WireNoDriver6e : public Layout
{
public:
    WireNoDriver6e (int N_p, double a = 1, double b = 3)
    {
        epc = epc6;

        double r_x = 0;
        double r_y = 0;
        
        for (int i=0; i<N_p; i++)
            addCell(r_x+i*(a+b), r_y, a);
     }
 };

bool epsilonEqual (double v, double w, double epsilon = 10E-8)
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
    QcaBond s1;
    s1.l = WireNoDriver2e(1);
    s1.constructBasis();
    BOOST_CHECK (s1.basis.size() == 6);

    QcaBond s2;
    s2.l = WireNoDriver2e(2);
    s2.constructBasis();
    BOOST_CHECK (s2.basis.size() == 36);

    QcaBond s3;
    s3.l = WireNoDriver2e(3);
    s3.constructBasis();
    BOOST_CHECK (s3.basis.size() == 216);
}

BOOST_AUTO_TEST_CASE ( test_qca_bond_system_for_some_parameters )
{
    /*
     * These values here are coming from tests/NotesAndTests.txt.
     */
    QcaBond s1;
    s1.l = WireNoDriver2e(1);
    s1.t = 0.2;
    s1.td = 0.04;
    s1.V0 = 10;
    //a = 1, b = 1
    s1.l = WireNoDriver2e(1, 1, 1);
    s1.constructBasis();
    s1.H.construct();
    s1.H.diagonalize();
    BOOST_CHECK (s1.energies().size() == 6);
    double eArray1[6] = {0.43, 0.43, 0.92, 1.08 , 1.28 , 1.28};
    std::vector<double> expected1(eArray1, eArray1+6);
    for (size_t i=0; i<expected1.size(); i++)
        BOOST_CHECK (epsilonEqual(expected1[i], s1.energies()(i), 0.01));

    // QcaBond s2;
    // s2.l = WireNoDriver2e(2);
    // s2.Vext = 0.1;
    // s2.t = 1;
    // s2.td = 0;
    // s2.V0 = 1000000;
    // //a = 0.001, b = 0.00175
    // s2.l = WireNoDriver2e(2, 0.001, 0.00175);
    // s2.update();
    // BOOST_CHECK (s2.energies().size() == 36);
    // double eArray2[5] = {2895.12, 2895.32, 2934.7, 2934.87, 2936.78};
    // std::vector<double> expected2(eArray2, eArray2+5);
    // for (size_t i=0; i<expected2.size(); i++)
    //     BOOST_CHECK (epsilonEqual(expected2[i], s2.energies()(i), 0.01));

    QcaBond s3;
    s3.l = Wire2e(1);
    double a = 1.0/250.0;
    double b = 1.75 * a;
    s3.t = 1;
    s3.td = 0.2;
    s3.V0 = 10.0 / a;
    //Pext = 0
    s3.l = Wire2e(1, a, b, 0);

    //Pext = 0.01
    s3.l = Wire2e(1, a, b, 0.01);
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000000, s3.P(0)), 0.29, 0.01));

    //Pext = 0.1
    s3.l = Wire2e(1, a, b, 0.1);
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000000, s3.P(0)), 0.93, 0.01));

    QcaBond s4;
    s4.l = Wire2e(3);
    s4.t = 1;
    s4.td = 0;
    s4.q = 0;
    s4.epsilonr = QCA_NATURAL_EPSILON_R; // natural units
    s4.lambdaD = 0; 
    s4.V0 = 1000;
    //a=0.01, b=0.023, Pext=1
    s4.l = Wire2e(3, 0.01, 0.023, 1);
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(0)), 0.670243, 1E-5));
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(1)), 0.462795, 1E-5));
    BOOST_CHECK (epsilonEqual(s4.measure(1, s4.P(2)), 0.295182, 1E-5));
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_fixed_charge_system )
{
    QcaFixedCharge s1;
    s1.l = WireNoDriver2e(1);
    s1.constructBasis();
    BOOST_CHECK (s1.basis.size() == 28);
    BOOST_CHECK (s1.basis.getRanges().size() == 3);

    QcaFixedCharge s2;
    s2.l = WireNoDriver2e(2);
    s2.constructBasis();
    BOOST_CHECK (s2.basis.size() == 28*28);
    BOOST_CHECK (s2.basis.getRanges().size() == 5);

#ifdef NDEBUG
    QcaFixedCharge s3;
    s3.l = WireNoDriver2e(3);
    s3.constructBasis();
    BOOST_CHECK (s3.basis.size() == 28*28*28);
    BOOST_CHECK (s3.basis.getRanges().size() == 7);
#endif

    QcaFixedCharge s4;
    s4.l = WireNoDriver6e(1);
    s4.constructBasis();
    BOOST_CHECK (s4.basis.size() == 28);

    QcaFixedCharge s5;
    s5.l = WireNoDriver6e(2);
    s5.constructBasis();
    BOOST_CHECK (s5.basis.size() == 28*28);
}

BOOST_AUTO_TEST_CASE ( test_qca_fixed_charge_system_for_some_parameters )
{
    double a,b,Pext;

    //see tests/NotesAndTests.txt
    QcaFixedCharge s1;
    s1.l = Wire2e(1);
    Pext = 0;
    s1.t = 1;
    s1.td = 0.2;
    a = 1.0/250.0;
    b = 1.75*a;
    s1.V0 = 10.0 / a;
    
    Pext = 0;
    s1.l = Wire2e(1, a, b, Pext);
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.N(0)), 2));
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.P(0)), 0));

    Pext = 0.1;
    s1.l = Wire2e(1, a, b, Pext);
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000, s1.P(0)), 0.895, 0.001));


    QcaFixedCharge s2;
    s2.l = Wire2e(2);
    Pext = 0;
    s2.t = 1;
    s2.td = 0.2;
    a = 1.0/250.0;
    b = 1.75*a;
    s2.V0 = 10.0 / a;
    
    Pext = 0;
    s2.l = Wire2e(2, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.N(0)), 2));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.N(1)), 2));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.P(0)), 0));
    BOOST_CHECK (epsilonEqual(s2.measure(1000, s2.P(1)), 0));

    Pext = 0.01;
    s2.l = Wire2e(2, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.662, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.013, 0.001));

    Pext = 0.2;
    s2.l = Wire2e(2, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.998, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.020, 0.001));

    Pext = 0;
    s2.t = 1;
    s2.td = 0;
    a = 1.0/100.0;
    b = 3*a;
    s2.V0 = 10.0 / a;
    
    Pext = 0.1;
    s2.l = Wire2e(2, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0.414, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), 0.379, 0.001));

    Pext = -0.8;
    s2.l = Wire2e(2, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), -0.944, 0.001));
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(1)), -0.850, 0.001));


    QcaFixedCharge s3;
    s3.l = Wire6e(1);
    Pext = 0;
    s3.t = 1;
    s3.td = 0.2;
    a = 1.0/250.0;
    b = 4*a;
    s3.V0 = 10.0 / a;
    
    Pext = 0;
    s3.l = Wire6e(1, a, b, Pext);
    s3.update();
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.N(0)), 6));
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.P(0)), 0));

    Pext = 0.1;
    s3.l = Wire6e(1, a, b, Pext);
    s3.update();
    BOOST_CHECK (s3.measure(1000, s3.P(0)) > 0.1);
    BOOST_CHECK (epsilonEqual(s3.measure(1000, s3.P(0)), 0.213, 0.001));


    QcaFixedCharge s4;
    s4.l = Wire6e(2);
    Pext = 0;
    s4.t = 1;
    s4.td = 0.2;
    a = 1.0/250.0;
    b = 5*a;
    s4.V0 = 10.0 / a;

    Pext = 0;
    s4.l = Wire6e(2, a, b, Pext);
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.N(0)), 6));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(0)), 0, 10E-8));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(1)), 0, 10E-8));

    Pext = 0.1;
    s4.l = Wire6e(2, a, b, Pext);
    s4.update();
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(0)), 0.268, 0.001));
    BOOST_CHECK (epsilonEqual(s4.measure(1000, s4.P(1)), 0.212, 0.001));
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_grandcanonical_system )
{
    QcaGrandCanonical s1;
    s1.l = WireNoDriver2e(1);
    s1.constructBasis();
    BOOST_CHECK (s1.basis.size() == 256);
    BOOST_CHECK (s1.basis.getRanges().size() == 1+2+3+4+5+4+3+2+1);

    QcaGrandCanonical s2;
    s2.l = WireNoDriver2e(2);
    s2.constructBasis();
    BOOST_CHECK (s2.basis.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( test_diagonalization_of_qca_grand_canonical_system )
{
    // TODO: right now this test does not really do anything

    QcaGrandCanonical s1;
    s1.l = WireNoDriver2e(1);
    s1.constructBasis();
    s1.H.construct();
    s1.H.diagonalize();
    BOOST_CHECK (s1.energies().size() == 256);

    QcaGrandCanonical s2;
    s2.l = WireNoDriver2e(2);
    s2.constructBasis();
    s2.H.construct();
    // this takes a _long_ time, but works
    // s2.H.diagonalize();
    // BOOST_CHECK (s2.H.eigenvalues.size() == 256*256);
}

BOOST_AUTO_TEST_CASE ( simple_sanity_checks_for_qca_grand_canonical_system )
{
    double a, b, Pext;

    QcaGrandCanonical s1;
    s1.l = Wire2e(1);
    Pext = 0;
    s1.t = 1;
    s1.td = 0;
    a = 1.0/160.0;
    b = 3*a;
    s1.V0 = 1000;
    s1.mu = 300;
    
    Pext = 0;
    s1.l = Wire2e(1, a, b, Pext);
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000000, s1.P(0)), 0));

    Pext = 0.1;
    s1.l = Wire2e(1, a, b, Pext);
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(1000000, s1.N(0)), 2));
    double P = s1.measure(100000, s1.P(0));
    BOOST_CHECK (P > 0.1 && P < 1);


    QcaGrandCanonical s2;
    s2.l = Wire6e(1);
    Pext = 0;
    s2.t = 1;
    s2.td = 0;
    a = 1.0/160.0;
    b = 4*a;
    s2.V0 = 1000;
    s2.mu = 1800;
    
    Pext = 0;
    s2.l = Wire6e(1, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.P(0)), 0));

    Pext = 1;
    s2.l = Wire6e(1, a, b, Pext);
    s2.update();
    BOOST_CHECK (epsilonEqual(s2.measure(1000000, s2.N(0)), 6));
    P = s2.measure(100000, s2.P(0));
    BOOST_CHECK (P > 0.1 && P < 1);
}

BOOST_AUTO_TEST_CASE ( performance_of_diagonalization_for_qca_fixed_charge )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;

    startCPUTime = std::clock();
    QcaFixedCharge s;
    s.l = Wire2e(2);
    s.constructBasis();
    s.H.construct();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction of two plaquet QCA fixed charge system: " 
              << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalize();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time to diagonalize a two plaquet QCA fixed charge system: " 
              << cpuTime << "s" << std::endl;
}

BOOST_AUTO_TEST_CASE ( test_scaling_of_parameters_for_grand_canonical_system )
{
    double a, b;

    QcaGrandCanonical s1;
    s1.l = WireNoDriver2e(1);
    s1.t = 1;
    s1.td = 0;
    a = 1.0/160.0;
    b = 3*a;
    s1.V0 = 1000;
    s1.mu = 300;
    s1.q = 0;

    s1.l = WireNoDriver2e(1, a, b);
    s1.update();
    double N1 = s1.measure(0.1, s1.N(0));
    double P1 = s1.measure(0.1, s1.P(0));

    QcaGrandCanonical s2;
    s2.l = WireNoDriver2e(1);
    s2.t = 1 * 10;
    s2.td = 0;
    a = 1.0/160.0 * 1.0/10.0;
    b = 3*a;
    s2.V0 = 1000 * 10;
    s2.mu = 300 * 10;
    s2.q = 0;

    s2.l = WireNoDriver2e(1, a, b);
    s2.update();
    double N2 = s2.measure(0.1 * 1.0/10.0, s2.N(0));
    double P2 = s2.measure(0.1 * 1.0/10.0, s2.P(0));

    BOOST_CHECK (epsilonEqual(N1, N2));
    BOOST_CHECK (epsilonEqual(P1, P2));


    // QcaGrandCanonical s3;
    // s3.l = WireNoDriver2e(1);
    // s3.Vext = 1;
    // s3.t = 1;
    // s3.td = 0;
    // a = 1.0/160.0;
    // b = 3*a;
    // s3.V0 = 1000;
    // s3.mu = 300;
    // s3.q = 0;

    // s3.l = WireNoDriver2e(1, a, b);
    // s3.update();
    // double N3 = s3.measure(0.1, s3.N(0));
    // double P3 = s3.measure(0.1, s3.P(0));

    // QcaGrandCanonical s4;
    // s4.l = WireNoDriver2e(1);
    // s4.Vext = 1 * 1000;
    // s4.t = 1 * 1000;
    // s4.td = 0;
    // a = 1.0/160.0 * 1.0/1000.0;
    // b = 3*a;
    // s4.V0 = 1000 * 1000;
    // s4.mu = 300 * 1000;
    // s4.q = 0;

    // s4.l = WireNoDriver2e(1, a, b);
    // s4.update();
    // double N4 = s4.measure(0.1 * 1.0/1000.0, s4.N(0));
    // double P4 = s4.measure(0.1 * 1.0/1000.0, s4.P(0));

    // BOOST_CHECK (epsilonEqual(N3, N4));
    // BOOST_CHECK (epsilonEqual(P3, P4));
}

BOOST_AUTO_TEST_CASE ( test_compensation_charge_2e_and_6e )
{
    double a, b, Pext;

    /*
     * Grand canonical, one plaquet: Should be four electrons per plaquet at 
     * mu = 0.
     */
    QcaGrandCanonical s1;
    s1.l = WireNoDriver2e(1);
    s1.t = 1;
    s1.td = 0;
    a = 1.0/160.0;
    b = 3*a;
    s1.V0 = 0;
    s1.mu = 0;
    s1.q = 1;

    s1.l = WireNoDriver2e(1, a, b);
    s1.update();
    BOOST_CHECK (epsilonEqual(s1.measure(10000, s1.N(0)), 4));

    /*
     * Fixed, two plaquets: For the case that doubly occupied states are
     * completely gapped out, the 2e per plaquet system without compensation
     * charges and the 6e per plaquet system with compensation charges should
     * behave identical.
     */
    QcaFixedCharge s2;
    s2.t = 1;
    s2.td = 0;
    a = 1.0/100.0;
    b = 2*a;
    s2.V0 = 10000;
    s2.mu = 0;
    s2.q = 0;

    s2.l = Wire2e(2, a, b, 1);
    s2.update();
    double P21 = s2.measure(10, s2.P(0));
    double P22 = s2.measure(10, s2.P(1));

    QcaFixedCharge s3;
    s3.t = 1;
    s3.td = 0;
    a = 1.0/100.0;
    b = 2*a;
    s3.V0 = 10000;
    s3.mu = 0;
    s3.q = 1;

    s3.l = Wire6e(2, a, b, 1);
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
    QcaGrandCanonical s4;
    s4.l = Wire6e(1);
    Pext = 1;
    s4.t = 1;
    s4.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s4.V0 = 1000;
    s4.mu = 1200;
    s4.q = 1;

    s4.l = Wire6e(1, a, b, Pext);
    s4.update();
    double P41 = s4.measure(10, s4.P(0));

    QcaGrandCanonical s5;
    s5.l = Wire2e(1);
    Pext = 1;
    s5.t = 1;
    s5.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s5.V0 = 1000;
    s5.mu = 200; //1000 less than s4.mu (because V0=1000)
    s5.q = 0;

    s5.l = Wire2e(1, a, b, Pext);
    s5.update();
    double P51 = s5.measure(10, s5.P(0));

    //std::cerr << "P41 = " << P41 << ", P51 = " << P51 << std::endl;

    BOOST_CHECK (epsilonEqual(P41, P51));

    /*
     * Fixed, two plaquets: The 2e and 6e dead plaquet systems should be
     * equivalent.
     */
    QcaFixedCharge s6;
    s6.l = Wire2e(2);
    Pext = 1;
    s6.t = 1;
    s6.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s6.V0 = 1000;
    s6.mu = 0;
    s6.q = 0;

    s6.l = Wire2e(2, a, b, Pext);
    s6.update();
    double P61 = s6.measure(10, s6.P(0));
    double P62 = s6.measure(10, s6.P(1));

    QcaFixedCharge s7;
    s7.l = Wire6e(2);
    Pext = 1;
    s7.t = 1;
    s7.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s7.V0 = 1000;
    s7.mu = 0;
    s7.q = 1;

    s7.l = Wire6e(2, a, b, Pext);
    s7.update();
    double P71 = s7.measure(10, s7.P(0));
    double P72 = s7.measure(10, s7.P(1));

    //std::cerr << "P61 = " << P61 << ", P62 = " << P62 << std::endl;
    //std::cerr << "P71 = " << P71 << ", P72 = " << P72 << std::endl;

    BOOST_CHECK (epsilonEqual(P61, P71));
    BOOST_CHECK (epsilonEqual(P62, P72));
}

BOOST_AUTO_TEST_CASE ( test_compensation_charge_2e_only )
{
    /*
     * Three tests related to compensations charges and two electrons per cell.
     * 
     * All tests use the bond model. The first two tests look at the expected
     * symmetries and inequalities of the occupancies of the dots without an
     * external potential (i.e. zero driving cell polarization), for a one cell
     * and two cell system, respectively. The third test asserts that setting a
     * finite input polarization yields expected results, namely that a net zero
     * cell charge (q = 1/2 for 2 electrons per cell) yields the best output
     * polarization.
     */

    // Test 1
    QcaBond s1;
    s1.l.nonuniformWire(1, 0.01, {0.01}, 0);
    s1.t = 1;
    s1.V0 = 1E6;
    
    std::vector<std::vector<std::vector<double>>> Ns;
    std::vector<double> qs = {0,0.25,0.5,0.75,1.0};
    for (double q: qs)
    {
        s1.q = q;
        s1.update();
        auto Ns_ = {s1.measureParticleNumber(0)};
        Ns.push_back(Ns_);
    }
    std::cerr << "Testing compensation charges for a one cell wire, "
              << "with zero input polarization." << std::endl;
    std::cerr << "   qs = ";
    for (auto q: qs) std::cerr << q << ", ";
    std::cerr << std::endl;
    for (auto Ns_: Ns)
    {
        std::cerr << "   One cell wire" << std::endl;
        for (auto N_: Ns_)
        {
            std::cerr << "      Ns = ";
            for (auto N: N_)
                std::cerr << N << " ";
            std::cerr << std::endl;
        }
    }
    for (auto Ns_: Ns)
    {
        BOOST_CHECK(epsilonEqual(Ns_[0][0], Ns_[0][1]));
        BOOST_CHECK(epsilonEqual(Ns_[0][2], Ns_[0][3]));
    }
    BOOST_CHECK(Ns[0][0][0] < 0.01);
    BOOST_CHECK(Ns[0][0][2] > 0.99);
    BOOST_CHECK(Ns[1][0][2] - Ns[1][0][0] > 0.05);
    BOOST_CHECK(epsilonEqual(Ns[2][0][0], 0.5));
    BOOST_CHECK(epsilonEqual(Ns[2][0][2], 0.5));
    BOOST_CHECK(Ns[3][0][0] - Ns[3][0][2] > 0.05);
    BOOST_CHECK(Ns[4][0][0] > 0.99);
    BOOST_CHECK(Ns[4][0][2] < 0.01);

    // Test 2
    s1.l.nonuniformWire(2, 0.01, {1000, 0.01}, 0);
    Ns = {};
    for (double q: qs)
    {
        s1.q = q;
        s1.update();
        auto Ns_ = {s1.measureParticleNumber(0),s1.measureParticleNumber(1)};
        Ns.push_back(Ns_);
    }
    std::cerr << "Testing compensation charges for a two cell wire, "
              << "with zero input polarization." << std::endl;
    std::cerr << "   qs = ";
    for (auto q: qs) std::cerr << q << ", ";
    std::cerr << std::endl;
    for (auto Ns_: Ns)
    {
        std::cerr << "   Two cell wire" << std::endl;
        for (auto N_: Ns_)
        {
            std::cerr << "      Ns = ";
            for (auto N: N_)
                std::cerr << N << " ";
            std::cerr << std::endl;
        }
    }
    for (auto Ns_: Ns)
    {
        BOOST_CHECK(epsilonEqual(Ns_[0][0], Ns_[0][1]));
        BOOST_CHECK(epsilonEqual(Ns_[0][2], Ns_[0][3]));
        BOOST_CHECK(epsilonEqual(Ns_[1][0], Ns_[1][1]));
        BOOST_CHECK(epsilonEqual(Ns_[1][2], Ns_[1][3]));
        BOOST_CHECK(epsilonEqual(Ns_[0][0], Ns_[1][2]));
        BOOST_CHECK(epsilonEqual(Ns_[0][2], Ns_[1][0]));
    }
    BOOST_CHECK(Ns[0][0][0] > 0.8);
    BOOST_CHECK(Ns[0][0][2] < 0.2);
    BOOST_CHECK(Ns[1][0][0] - Ns[1][0][2] > 0.001);
    BOOST_CHECK(epsilonEqual(Ns[2][0][0], 0.5));
    BOOST_CHECK(epsilonEqual(Ns[2][0][2], 0.5));
    BOOST_CHECK(Ns[3][0][2] - Ns[3][0][0] > 0.001);
    BOOST_CHECK(Ns[4][0][2] > 0.8);
    BOOST_CHECK(Ns[4][0][0] < 0.2);

    // Test 3
    s1.l.nonuniformWire(2, 0.02, {0.04, 0.04}, 1);
    Ns = {};
    std::vector<std::vector<double>> Ps;
    for (double q: qs)
    {
        s1.q = q;
        s1.update();
        auto Ns_ = {s1.measureParticleNumber(0),s1.measureParticleNumber(1)};
        auto Ps_ = {s1.measurePolarization(0),s1.measurePolarization(1)};
        Ns.push_back(Ns_);
        Ps.push_back(Ps_);
    }
    std::cerr << "Testing compensation charges for a two cell wire, "
              << "with finite input polarization." << std::endl;
    std::cerr << "   qs = ";
    for (auto q: qs) std::cerr << q << ", ";
    std::cerr << std::endl;
    // for (auto Ns_: Ns)
    // {
    //     std::cerr << "   Two cell wire" << std::endl;
    //     for (auto N_: Ns_)
    //     {
    //         std::cerr << "      Ns = ";
    //         for (auto N: N_)
    //             std::cerr << N << " ";
    //         std::cerr << std::endl;
    //     }
    // }
    for (auto Ps_: Ps)
    {
        std::cerr << "   Two cell wire" << std::endl;
        for (auto P: Ps_)
            std::cerr << "      P = " << P << std::endl;
    }
    BOOST_CHECK(epsilonEqual(Ps[0][0], Ps[4][0]));
    BOOST_CHECK(epsilonEqual(Ps[0][1], Ps[4][1]));
    BOOST_CHECK(epsilonEqual(Ps[1][0], Ps[3][0]));
    BOOST_CHECK(epsilonEqual(Ps[1][1], Ps[3][1]));
    BOOST_CHECK(Ps[2][0] > Ps[1][0]);
    BOOST_CHECK(Ps[2][1] > Ps[1][1]);
    BOOST_CHECK(Ps[1][0] > Ps[0][0]);
    BOOST_CHECK(Ps[1][1] > Ps[0][1]);
}

BOOST_AUTO_TEST_CASE ( compare_polarization_and_polarization2 )
{
    double a, b, Pext;

    QcaGrandCanonical s1;
    s1.l = Wire6e(1);
    Pext = 1;
    s1.t = 1;
    s1.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s1.V0 = 1000;
    s1.mu = 1200;
    s1.q = 1;

    s1.l = Wire6e(1, a, b, Pext);
    s1.update();
    double P11 = s1.measure(10, s1.P(0));
    double P12 = s1.measurePolarization2(10, 0);

    //std::cerr << "P11 = " << P11 << ", P12 = " << P12 << std::endl;
    BOOST_CHECK (epsilonEqual(s1.measure(10, s1.N(0)), 6));
    BOOST_CHECK (epsilonEqual(P11, P12, 1e-5));

    QcaGrandCanonical s2;
    s2.l = Wire6e(1);
    Pext = 1;
    s2.t = 1;
    s2.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s2.V0 = 10;
    s2.mu = 200;
    s2.q = 1;

    s2.l = Wire6e(1, a, b, Pext);
    s2.update();
    double P21 = s2.measure(10, s2.P(0));
    double P22 = s2.measurePolarization2(10, 0);

    //std::cerr << "P21 = " << P21 << ", P22 = " << P22 << std::endl;
    BOOST_CHECK (epsilonEqual(s2.measure(10, s2.N(0)), 6));
    BOOST_CHECK (P21 > 0.9);
    BOOST_CHECK (P22 < 0.1);

    QcaFixedCharge s3;
    s3.l = Wire6e(2);
    Pext = 1;
    s3.t = 1;
    s3.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s3.V0 = 1000;
    s3.mu = 0;
    s3.q = 1;

    s3.l = Wire6e(2, a, b, Pext);
    s3.update();
    double P311 = s3.measure(10, s3.P(0));
    double P312 = s3.measure(10, s3.P(1));
    double P321= s3.measurePolarization2(10, 0);
    double P322= s3.measurePolarization2(10, 1);

    BOOST_CHECK (epsilonEqual(P311, P321, 1e-5));
    BOOST_CHECK (epsilonEqual(P312, P322, 1e-5));

    QcaFixedCharge s4;
    s4.l = Wire6e(2);
    Pext = 1;
    s4.t = 1;
    s4.td = 0;
    a = 1.0/100.0;
    b = 4*a;
    s4.V0 = 10;
    s4.mu = 0;
    s4.q = 1;

    s4.l = Wire6e(2, a, b, Pext);
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
    l1.addSite(0,0);
    l1.addSite(0,1);
    l1.addSite(1,1);
    l1.addSite(1,0);
    l1.addCharge(-1,1,1);
    BOOST_CHECK (l1.r(0,1) == 1);
    BOOST_CHECK (l1.r(1,3) == std::sqrt(2));
    BOOST_CHECK (l1.r_charge_dot(0,3) == std::sqrt(5));

    Layout l2;
    l2.addCell(0,0,0.2);
    l2.addCell(0,2,0.2);
    l2.addDriverCell(1,0,0.2,1);
    l2.addSite(1,1);
    l2.addCharge(0,-1,0.5);
    BOOST_CHECK(epsilonEqual(l2.r(0,3), 0.2));
    BOOST_CHECK(epsilonEqual(l2.r(4,6), std::sqrt(0.2*0.2+0.2*0.2)));
    BOOST_CHECK(epsilonEqual(l2.r(0,4), 2));
    BOOST_CHECK(epsilonEqual(l2.r(1,7), std::sqrt(1.8*1.8+0.2*0.2)));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(0,0), 1));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(2,1), 1.2));
    BOOST_CHECK(epsilonEqual(l2.r(0,8), std::sqrt(2)));
    BOOST_CHECK(epsilonEqual(l2.r(6,8), std::sqrt(1.2*1.2+0.8*0.8)));
    BOOST_CHECK(epsilonEqual(l2.r_charge_dot(4,3), std::sqrt(1+0.2*0.2)));
    BOOST_CHECK(l2.N_sites() == 9);
    BOOST_CHECK(l2.N_charges() == 5);

    Layout l3;
    l3.wire(4, 0.2, 0.4, 0.7);
    BOOST_CHECK(l3.N_sites() == 16);
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
    l4.nonuniformWire(4, 0.2, bs, -0.2, epc6);
    BOOST_CHECK(l4.N_sites() == 16);
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
    l5.nonuniformWire(3, 0.1, bs, 0);
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
    l6.addDriverCell(0,0, 1, 0.3, epc6);
    BOOST_CHECK(epsilonEqual(
                0.5*(l6.charge(0)+l6.charge(2)-l6.charge(1)-l6.charge(3)), 
                0.3
    ));

    Layout l7;
    l7.addDriverCell(0,0, 1, -0.12, epc6);
    BOOST_CHECK(epsilonEqual(
                0.5*(l7.charge(0)+l7.charge(2)-l7.charge(1)-l7.charge(3)), 
                -0.12
    ));
}

BOOST_AUTO_TEST_CASE ( reuse_same_system_multiple_times_with_different_layouts )
{
    // This test case also documents the new preferred usage of the QCA classes.
    // note: default is natural / dimensionless units
    QcaBond s1;
    s1.l = Wire2e(2, 0.01, 0.02, 1);
    s1.beta = 1;
    s1.update();
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                0.863919,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(1),
                0.603095,
                1E-5
                ));
    
    // for changed beta, s1.update() is not necessary
    s1.beta = 2;
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                0.986214,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(1),
                0.813721,
                1E-5
                ));

    // change layout, for example, inter-cell spacing
    s1.l = Wire2e(2, 0.02, 0.02, 1);
    s1.update();
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                0.978439,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(1),
                0.00399499,
                1E-5
                ));

    // change layout completely: different number of cells
    s1.l = Wire2e(3, 0.03, 0.06, 1);
    s1.update();
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                0.668657,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(1),
                0.455502,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(2),
                0.115793,
                1E-5
                ));
    
    // again, change layout completely
    s1.l.wire(1, 0.01, 0.03, 0.4);
    s1.update();
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                0.234819,
                1E-5
                ));

    // test fixed charge system as well
    QcaFixedCharge s2;
    s2.l = Wire2e(2, 0.01, 0.02, 0.1);
    s2.V0 = 1000;
    s2.beta = 100;
    s2.update();
    BOOST_CHECK (epsilonEqual(
                s2.measurePolarization(0),
                0.940904,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s2.measurePolarization(1),
                0.765352,
                1E-5
                ));
    
    // chage to 6 electrons per cell instead of two, and also different number
    // of cells
    s2.l.epc = epc6;
    s2.l.wire(1, 0.01, 0.04, 0.1, epc6);
    s2.update();
    BOOST_CHECK (epsilonEqual(
                s2.measurePolarization(0),
                0.0347185,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s2.measure(s2.beta, s2.N(0)),
                6,
                1E-5
                ));
    
    // now two cells instead of one cell
    s2.l = Wire6e(2, 0.01, 0.04, 0.1);
    s2.update();
    BOOST_CHECK (epsilonEqual(
                s2.measurePolarization(0),
                0.0718588,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s2.measurePolarization(1),
                0.00265602,
                1E-5
                ));
}

BOOST_AUTO_TEST_CASE ( limits_of_nonuniform_wire )
{
    // A non-uniform with uniform inter-cell spacing is just a uniform wire
    std::vector<double> bs1;
    bs1.push_back(0.02);
    bs1.push_back(0.02);
    bs1.push_back(0.02);
    
    QcaBond s1;
    s1.l.nonuniformWire(3, 0.01, bs1, 1);
    s1.beta = 1;
    s1.update();
    
    QcaBond s2;
    s2.l.wire(3, 0.01, 0.02, 1);
    s2.beta = 1;
    s2.update();

    //std::cerr 
    //    << s1.measurePolarization(0) << "  " 
    //    << s1.measurePolarization(1) << "  " 
    //    << s1.measurePolarization(2) << std::endl;
   
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(0),
                s2.measurePolarization(0),
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(1),
                s2.measurePolarization(1),
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s1.measurePolarization(2),
                s2.measurePolarization(2),
                1E-5
                ));

    // Same, but for fixed charge and 6 electrons per cell
    std::vector<double> bs2;
    bs2.push_back(0.035);
    bs2.push_back(0.035);
    
    QcaFixedCharge s3;
    s3.l.epc = epc6;
    s3.l.nonuniformWire(2, 0.01, bs2, 1, epc6);
    s3.V0 = 1000;
    s3.beta = 1;
    s3.update();

    QcaFixedCharge s4;
    s4.l.epc = epc6;
    s4.l.wire(2, 0.01, 0.035, 1, epc6);
    s4.V0 = 1000;
    s4.beta = 1;
    s4.update();

    //std::cerr 
    //    << s3.measurePolarization(0) << "  " 
    //    << s3.measurePolarization(1) << std::endl;

    BOOST_CHECK (epsilonEqual(
                s3.measurePolarization(0),
                s4.measurePolarization(0),
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s3.measurePolarization(1),
                s4.measurePolarization(1),
                1E-5
                ));

    // In the limit of a very large inter-cell spacing, we should see no
    // polarization response.
    // First we put the driver cell far away.
    std::vector<double> bs3;
    bs3.push_back(10);
    bs3.push_back(0.02);
    bs3.push_back(0.02);

    QcaBond s5;
    s5.l.nonuniformWire(3, 0.01, bs3, 1);
    s5.beta = 1;
    s5.update();

    BOOST_CHECK (epsilonEqual(
                s5.measurePolarization(0),
                0,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s5.measurePolarization(1),
                0,
                1E-5
                ));
    BOOST_CHECK (epsilonEqual(
                s5.measurePolarization(2),
                0,
                1E-5
                ));

    // Now we move the right-most cell (the output cell) far away.
    std::vector<double> bs4;
    bs4.push_back(0.02);
    bs4.push_back(0.02);
    bs4.push_back(10);

    QcaBond s6;
    s6.l.nonuniformWire(3, 0.01, bs4, 1);
    s6.beta = 1;
    s6.update();

    BOOST_CHECK ( s6.measurePolarization(0) > 0.3);
    BOOST_CHECK ( s6.measurePolarization(1) > 0.3);
    BOOST_CHECK (epsilonEqual(
                s6.measurePolarization(2),
                0,
                1E-5
                ));
}
