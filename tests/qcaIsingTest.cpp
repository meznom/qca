#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE qca test
#include <boost/test/unit_test.hpp>
#include "qca.hpp"

bool equal(double a, double b, double epsilon=1E-8)
{
    return std::abs(a-b) < epsilon;
}

BOOST_AUTO_TEST_CASE ( test_construct_qca_ising_system )
{
    QcaIsing s;

    s.l.wire(1,1,1,1);
    s.update();
    BOOST_CHECK (s.basis.size() == 2);

    s.l.wire(2,1,1,1);
    s.update();
    BOOST_CHECK (s.basis.size() == 4);

    s.l.wire(6,1,1,1);
    s.update();
    BOOST_CHECK (s.basis.size() == 64);
}

BOOST_AUTO_TEST_CASE ( test_sigma_operator )
{
    const size_t N_c = 5;
    
    QcaIsing s;
    s.l.wire(N_c,1,1,1);
    s.update();

    State state("10010");
    std::vector<double> expected_spins = {1,-1,-1,1,-1};

    SparseVector<double> a(s.basis.size());
    a.reserve(1);
    a.insert(s.basis(state)) = 1;

    std::vector<double> spins(N_c);
    for (size_t j=0; j<N_c; j++)
    {
        SMatrix m = a.adjoint() * s.sigma(j) * a;
        assert(m.size() == 1);
        spins[j] = m.coeffRef(0,0);
        // std::cerr << "Spin " << j+1 << ": " << spins[j] << std::endl;
    }
    BOOST_CHECK (expected_spins == spins);
}

BOOST_AUTO_TEST_CASE ( test_spin_and_polarization_are_the_same )
{
    QcaIsing s;
    s.beta = 1;
    s.t = 1E-3;
    s.l.wire(3,0.01,0.02,1);
    s.update();

    std::vector<double> spins = {s.measureSpin(0), 
                                 s.measureSpin(1), 
                                 s.measureSpin(2)};
    std::vector<double> Ps = {s.measurePolarization(0), 
                              s.measurePolarization(1), 
                              s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        // std::cerr << "Spin[" << i << "] = " << spins[i] 
        //           << ", Polarization[" << i << "] = " << Ps[i] << std::endl;
        BOOST_CHECK (equal(spins[i], Ps[i]));
    }
}

BOOST_AUTO_TEST_CASE ( test_simple_physical_limits )
{
    QcaIsing s;
    s.beta = 1;
    s.t = 1E-3;
    s.l.wire(3,0.01,0.02,0);
    s.update();

    /*
     * Different input polarizations.
     */
    // P^D = 0
    std::vector<double> Ps1 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
        BOOST_CHECK (equal(0, Ps1[i]));
    // P^D = 1
    s.l.wire(3,0.01,0.02,1);
    s.update();
    std::vector<double> Ps2 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
        BOOST_CHECK (Ps2[i]>0.5);
    for (size_t i=1; i<3; i++)
        BOOST_CHECK (Ps2[i]<Ps2[i-1]);
    // P^D = -1
    s.l.wire(3,0.01,0.02,-1);
    s.update();
    std::vector<double> Ps3 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        BOOST_CHECK (Ps3[i]<-0.5);
        BOOST_CHECK (equal(Ps2[i], -Ps3[i]));
    }
    for (size_t i=1; i<3; i++)
        BOOST_CHECK (Ps3[i]>Ps3[i-1]);

    /*
     * Different temperatures.
     */
    s.l.wire(3,0.01,0.02,1);
    s.update();
    // High temperature
    s.beta = 0.01;
    std::vector<double> Ps4 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    // Medium temperature
    s.beta = 1;
    std::vector<double> Ps5 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    // Low temperature
    s.beta = 100;
    std::vector<double> Ps6 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        BOOST_CHECK (Ps4[i] < Ps5[i]);
        BOOST_CHECK (Ps5[i] < Ps6[i]);
        BOOST_CHECK (Ps4[i] < 0.1);
    }

    /*
     * Hopping versus Coulomb potential.
     *
     * Large hopping should wash out any polarization.
     */
    // Small hopping
    s.l.wire(3,0.01,0.02,1);
    s.t = 1E-3;
    s.beta = 1;
    s.update();
    std::vector<double> Ps7 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    // Medium hopping
    s.t = 1;
    s.update();
    std::vector<double> Ps8 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    // Large hopping, on the order of V_1
    s.t = 100;
    s.update();
    std::vector<double> Ps9 = {s.measurePolarization(0), 
                               s.measurePolarization(1), 
                               s.measurePolarization(2)};
    // Medium hopping, small V_1
    s.l.wire(3,0.1,0.2,1);
    s.t = 1;
    s.update();
    std::vector<double> Ps10 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    // Medium hopping, medium V_1
    s.l.wire(3,0.01,0.02,1);
    s.t = 1;
    s.update();
    std::vector<double> Ps11 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    // Medium hopping, large V_1
    s.l.wire(3,0.001,0.002,1);
    s.t = 1;
    s.update();
    std::vector<double> Ps12 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        BOOST_CHECK (Ps7[i] > Ps8[i]);
        BOOST_CHECK (Ps8[i] > Ps9[i]);
        BOOST_CHECK (Ps9[i] < 1E-2);
        BOOST_CHECK (Ps10[i] < Ps11[i]);
        BOOST_CHECK (Ps11[i] < Ps12[i]);
        BOOST_CHECK (Ps12[i] > 0.95);
    }

    /*
     * Larger inter-cell spacing means smaller polarization.
     */
    s.l.wire(3,0.01,0.02,1);
    s.t = 1E-3;
    s.beta = 1;
    s.update();
    std::vector<double> Ps13 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    s.l.wire(3,0.01,0.04,1);
    s.update();
    std::vector<double> Ps14 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    s.l.wire(3,0.01,1,1);
    s.update();
    std::vector<double> Ps15 = {s.measurePolarization(0), 
                                s.measurePolarization(1), 
                                s.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        BOOST_CHECK (Ps13[i] > Ps14[i]);
        BOOST_CHECK (Ps14[i] > Ps15[i]);
        BOOST_CHECK (Ps15[i] < 1E-3);
    }
}

BOOST_AUTO_TEST_CASE ( test_tprime_gets_set_and_calculated_correctly )
{
    /*
     * t^{\prime} = \frac{ 4 t^2 }{ \Delta V }
     *            = \frac{ 8 t^2 }{ (2 - \sqrt{2}) V_1 }
     *            = \frac{ 4 t^2 }{ 0.29289 V_1 }
     */
    double t;
    double V1;

    QcaIsing s;
    s.q = 0.5;
    s.beta = 1;

    t = 1;
    V1 = 100;
    
    s.t = t;
    s.l.wire(3,1.0/V1,1,1);
    s.update();
    
    // std::cerr << "t' = " << s.tprime << std::endl;
    BOOST_CHECK(equal(s.tprime, 0.13656854249492380195));

    t = 3;
    V1 = 40;
    s.t = t;
    s.l.wire(3,1.0/V1,1,1);
    s.update();
    
    // std::cerr << "t' = " << s.tprime << std::endl;
    BOOST_CHECK(equal(s.tprime, 3.07279220613578554391));

    t = 10;
    V1 = 200;
    s.t = t;
    s.l.wire(3,1.0/V1,1,1);
    s.update();
    
    // std::cerr << "t' = " << s.tprime << std::endl;
    BOOST_CHECK(equal(s.tprime, 6.82842712474619009758));
}

BOOST_AUTO_TEST_CASE ( test_bond_and_ising_models_are_the_same_in_some_limits )
{
    QcaIsing s_i;
    s_i.q = 0.5;
    s_i.beta = 1;
    s_i.t = 1;
    s_i.l.wire(3,0.001,0.005,1);
    s_i.update();

    QcaBond s_b;
    s_b.q = 0.5;
    s_b.beta = 1;
    s_b.t = 1;
    s_b.l.wire(3,0.001,0.005,1);
    s_b.update();

    std::vector<double> Ps_i = {s_i.measurePolarization(0), 
                                s_i.measurePolarization(1), 
                                s_i.measurePolarization(2)};
    std::vector<double> Ps_b = {s_b.measurePolarization(0), 
                                s_b.measurePolarization(1), 
                                s_b.measurePolarization(2)};
    for (size_t i=0; i<3; i++)
    {
        std::cerr << "Ps_i[" << i << "] = " << Ps_i[i] 
                  << ", Ps_b[" << i << "] = " << Ps_b[i] << std::endl;
    }
    // TODO
}
