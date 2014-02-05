#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE system test
#include <boost/test/unit_test.hpp>

#include <ctime>
#include "basis.hpp"
#include "system.hpp"

template<class System>
class HubbardHamiltonian : public Hamiltonian<System>
{
public:
    HubbardHamiltonian (System& s_) 
        : Hamiltonian<System>(s_), H(Hamiltonian<System>::H), s(s_)
    {}

    void construct ()
    {
        H = SMatrix(s.basis.size(), s.basis.size());
        H.setZero();
        for (size_t i=0; i<s.N_sites; i++)
        {
            for (int spin=0; spin<2; spin++)
            {
                H += s.t * s.creator(2*i+spin) * s.annihilator(2* ((i+1)%s.N_sites) +spin);
                H += s.t * s.creator(2* ((i+1)%s.N_sites) +spin) * s.annihilator(2*i+spin);
            }
            H += s.U * s.creator(2*i) * s.annihilator(2*i) *
                       s.creator(2*i+1) * s.annihilator(2*i+1);
        }
    }

    SMatrix& H;
    System& s;
};

typedef BasicSystem BaseSystem;
class HubbardSystem : public BaseSystem
{
public:
    HubbardSystem (size_t N_sites_, bool exploitSymmetries_ = false) 
        : BaseSystem(2*N_sites_), N_sites(N_sites_), H(*this), 
          exploitSymmetries(exploitSymmetries_), measure(*this), measureBS(*this),
          t(1), U(10)
    {
        if (exploitSymmetries)
        {
            //TODO: technically there's a problem here, because N and S get
            //destroyed before basis gets destroyed
            basis.addSymmetryOperator(&N);
            basis.addSymmetryOperator(&S);
        }
        BaseSystem::construct();
    }

    void update()
    {
        H.construct();
        H.diagonalize();
    }

    size_t N_sites;
    HubbardHamiltonian<HubbardSystem> H;
    bool exploitSymmetries;
    ParticleNumberSymmetryOperator N;
    SpinSymmetryOperator S;
    EnsembleAverage<HubbardSystem> measure;
    EnsembleAverageBySectors<HubbardSystem> measureBS;
    double t, U;
};

template<class System>
class DoubleOccupancy
{
public:
    DoubleOccupancy (const System& s_)
    : s(s_)
    {}

    SMatrix operator() (size_t site) const
    {
        return 
            s.creator(2*site) * s.annihilator(2*site) *
            s.creator(2*site+1) * s.annihilator(2*site+1);
    }
private:
    const System& s;
};

bool operator== (SMatrix& sm1, SMatrix& sm2)
{
    if (sm1.rows() != sm2.rows() || sm1.cols() != sm2.cols())
        return false;
    for (int j=0; j<sm1.cols(); j++)
        for (int i=0; i<sm1.rows(); i++)
            if (sm1.coeffRef(i,j) != sm2.coeffRef(i,j))
                return false;
    return true;
}

class EpsilonEqualPred
{
public:
    EpsilonEqualPred (double v_, double epsilon_ = 10E-10)
        : v(v_), epsilon(epsilon_)
    {}

    bool operator() (double w) const
    {
        return fabs(v - w) < epsilon;
    }

private:
    double v;
    double epsilon;
};

bool epsilonEqual (double v, double w, double epsilon = 10E-10)
{
    return fabs(v-w) < epsilon;
}

BOOST_AUTO_TEST_CASE ( construct_system_without_symmetries )
{
    HubbardSystem s(2);
    s.H.construct();
    BOOST_CHECK (s.basis.getRanges().size() == 1);
    
    s.H.diagonalizeNoSymmetries();
    DVector eigenvalues1 = s.H.eigenvalues();
    SMatrix eigenvectors1 = s.H.eigenvectors();
    
    s.H.diagonalizeUsingSymmetries();
    DVector eigenvalues2 = s.H.eigenvalues();
    SMatrix eigenvectors2 = s.H.eigenvectors();

    s.H.diagonalizeUsingSymmetriesBySectors();
    DVector eigenvalues3 = s.H.eigenvalues();
    SMatrix eigenvectors3 = s.H.eigenvectors();

    BOOST_CHECK (eigenvalues1 == eigenvalues2);
    BOOST_CHECK (eigenvectors1.size() == eigenvectors2.size());
    BOOST_CHECK (eigenvalues1 == eigenvalues3);
    BOOST_CHECK (eigenvectors1.size() == eigenvectors3.size());
}

BOOST_AUTO_TEST_CASE ( construct_system_with_symmetries )
{
    HubbardSystem s(4, true);
    s.H.construct();
    BOOST_CHECK (s.basis.getRanges().size() > 1);

    s.H.diagonalizeNoSymmetries();
    DVector eigenvalues1 = s.H.eigenvalues();
    auto& evs1 = s.H.eigenvectors();
    DMatrix eigenvectors1 = evs1.block(0,0,evs1.cols(),evs1.cols());
    std::vector<double> ev1;
    for (int i=0; i<eigenvalues1.size(); i++)
        ev1.push_back(eigenvalues1(i));
    
    s.H.diagonalizeUsingSymmetries();
    DVector eigenvalues2 = s.H.eigenvalues();
    auto& evs2 = s.H.eigenvectors();
    DMatrix eigenvectors2 = evs2.block(0,0,evs2.cols(),evs2.cols());
    
    std::vector<double> ev2;
    for (int i=0; i<eigenvalues2.size(); i++)
        ev2.push_back(eigenvalues2(i));

    s.H.diagonalizeUsingSymmetriesBySectors();
    DVector eigenvalues3 = s.H.eigenvalues();
    auto& evs3 = s.H.eigenvectors();
    DMatrix eigenvectors3 = evs3.block(0,0,evs3.cols(),evs3.cols());
    
    std::vector<double> ev3;
    for (int i=0; i<eigenvalues3.size(); i++)
        ev3.push_back(eigenvalues3(i));

    BOOST_CHECK (ev1.size() == ev2.size());
    BOOST_CHECK (ev1.size() == ev3.size());

    //We only check if both methods yield the same eigenvalues. Checking
    //whether the eigenvectors are equivalent is too complicated.
    for (size_t i=0; i<ev1.size(); i++)
    {
        BOOST_CHECK (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i])) != ev2.end());
        BOOST_CHECK (std::find_if(ev3.begin(), ev3.end(), EpsilonEqualPred(ev1[i])) != ev3.end());
    }
}

BOOST_AUTO_TEST_CASE ( construct_and_diagonalize_system_multiple_times )
{
    HubbardSystem s(4, true);
    s.H.construct();
    size_t n_b1 = s.basis.size();
    size_t n_r1 = s.basis.getRanges().size();
    
    s.H.diagonalizeUsingSymmetries();
    DVector ev1 = s.H.eigenvalues();

    s.H.construct();
    size_t n_b2 = s.basis.size();
    size_t n_r2 = s.basis.getRanges().size();
    BOOST_CHECK (n_b1 == n_b2);
    BOOST_CHECK (n_r1 == n_r2);
   
    s.H.diagonalizeUsingSymmetries();
    DVector ev2 = s.H.eigenvalues();

    BOOST_CHECK (ev1 == ev2);
}

BOOST_AUTO_TEST_CASE ( test_BySectors_eigenvalue_and_eigenvectors_accessors )
{
    HubbardSystem s(4, true);
    s.H.construct();
   

    s.H.diagonalizeNoSymmetries();
    std::vector<DVector> eigenvalues1 = s.H.eigenvaluesBySectors();

    BOOST_CHECK (eigenvalues1.size() == 1);
    BOOST_CHECK (eigenvalues1[0].size() > 1);

    std::vector<double> ev1;
    for (int i=0; i<eigenvalues1[0].size(); i++)
        ev1.push_back(eigenvalues1[0](i));


    s.H.diagonalizeUsingSymmetries();
    std::vector<DVector> eigenvalues2 = s.H.eigenvaluesBySectors();
    
    BOOST_CHECK (eigenvalues2.size() == 1);
    BOOST_CHECK (eigenvalues2[0].size() > 1);

    std::vector<double> ev2;
    for (int i=0; i<eigenvalues2[0].size(); i++)
        ev2.push_back(eigenvalues2[0](i));


    s.H.diagonalizeUsingSymmetriesBySectors();
    std::vector<DVector> eigenvalues3 = s.H.eigenvaluesBySectors();

    BOOST_CHECK (eigenvalues3.size() > 1);
    BOOST_CHECK (eigenvalues3[0].size() > 0);

    std::vector<double> ev3;
    for (size_t j=0; j<eigenvalues3.size(); j++)
        for (int i=0; i<eigenvalues3[j].size(); i++)
            ev3.push_back(eigenvalues3[j](i));


    BOOST_CHECK (ev1.size() == ev2.size());
    BOOST_CHECK (ev1.size() == ev3.size());

    for (size_t i=0; i<ev1.size(); i++)
    {
        BOOST_CHECK (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i])) != ev2.end());
        BOOST_CHECK (std::find_if(ev3.begin(), ev3.end(), EpsilonEqualPred(ev1[i])) != ev3.end());
    }


    // test that eigenvectors get properly reset and reconstructed
    s.H.diagonalizeNoSymmetries();
    std::vector<DVector> eigenvalues4 = s.H.eigenvaluesBySectors();

    BOOST_CHECK (eigenvalues4.size() == 1);
    BOOST_CHECK (eigenvalues4[0].size() > 1);

    std::vector<double> ev4;
    for (int i=0; i<eigenvalues4[0].size(); i++)
        ev4.push_back(eigenvalues4[0](i));

    s.t = 1000;
    s.U = 20000;
    s.H.construct();
    s.H.diagonalizeUsingSymmetries();

    std::vector<DVector> eigenvalues5 = s.H.eigenvaluesBySectors();

    BOOST_CHECK (eigenvalues5.size() == 1);
    BOOST_CHECK (eigenvalues5[0].size() > 1);
    
    std::vector<double> ev5;
    for (int i=0; i<eigenvalues5[0].size(); i++)
        ev5.push_back(eigenvalues5[0](i));

    bool in = false;
    for (size_t i=0; i<ev4.size(); i++)
        if (ev4[i] > 10E-10 && 
            std::find_if(ev5.begin(), ev5.end(), EpsilonEqualPred(ev4[i])) != ev5.end())
        {
            in = true;
        }
    BOOST_CHECK (in == false);
}

BOOST_AUTO_TEST_CASE ( measure_double_occupancy )
{
    HubbardSystem s(4, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.H.construct();
    s.H.diagonalizeNoSymmetries();
    double doLowT1 = s.measure(1000, DO(0));
    double doMidT1 = s.measure(1, DO(0));
    double doHighT1 = s.measure(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetries();
    double doLowT2 = s.measure(1000, DO(0));
    double doMidT2 = s.measure(1, DO(0));
    double doHighT2 = s.measure(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetriesBySectors();
    double doLowT3 = s.measure(1000, DO(0));
    double doMidT3 = s.measure(1, DO(0));
    double doHighT3 = s.measure(0.00001, DO(0));

    BOOST_CHECK (epsilonEqual(doLowT1, 0, 10E-2));
    BOOST_CHECK (epsilonEqual(doHighT1, 0.25, 10E-2));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT2));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT2));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT2));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT3));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT3));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT3));
}

BOOST_AUTO_TEST_CASE ( test_ensembleAverageBySectors )
{
    HubbardSystem s(4, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.H.construct();
    s.H.diagonalizeNoSymmetries();
    double doLowT1 = s.measureBS(1000, DO(0));
    double doMidT1 = s.measureBS(1, DO(0));
    double doHighT1 = s.measureBS(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetries();
    double doLowT2 = s.measureBS(1000, DO(0));
    double doMidT2 = s.measureBS(1, DO(0));
    double doHighT2 = s.measureBS(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetriesBySectors();
    double doLowT3 = s.measureBS(1000, DO(0));
    double doMidT3 = s.measureBS(1, DO(0));
    double doHighT3 = s.measureBS(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetriesBySectors();
    double doLowT4 = s.measure(1000, DO(0));
    double doMidT4 = s.measure(1, DO(0));
    double doHighT4 = s.measure(0.00001, DO(0));

    s.H.diagonalizeUsingSymmetries();
    double doLowT5 = s.measure(1000, DO(0));
    double doMidT5 = s.measure(1, DO(0));
    double doHighT5 = s.measure(0.00001, DO(0));

    BOOST_CHECK (epsilonEqual(doLowT1, 0, 10E-2));
    BOOST_CHECK (epsilonEqual(doHighT1, 0.25, 10E-2));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT2));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT2));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT2));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT3));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT3));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT3));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT4));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT4));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT4));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT5));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT5));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT5));
}

BOOST_AUTO_TEST_CASE ( compare_performance_diagonalize_and_diagonalizeUsingSymmetries_and_diagonalizeUsingSymmetriesBySectors )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    size_t sites;
#ifdef NDEBUG
    sites = 5;
#else
    sites = 4;
#endif

    startCPUTime = std::clock();
    HubbardSystem s(sites, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.H.construct();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeNoSymmetries();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalize of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetries();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetries of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetriesBySectors();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetriesBySectors of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;
}

BOOST_AUTO_TEST_CASE ( compare_performance_diagonalizeUsingSymmetries_and_diagonalizeUsingSymmetriesBySectors )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    size_t sites;
#ifdef NDEBUG
    sites = 7;
#else
    sites = 5;
#endif

    startCPUTime = std::clock();
    HubbardSystem s(sites, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.H.construct();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for construction of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetries();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetries of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.H.diagonalizeUsingSymmetriesBySectors();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeUsingSymmetriesBySectors of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;
}

BOOST_AUTO_TEST_CASE ( compare_performance_ensembleAverage_and_ensembleAverageBySectors )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    size_t sites;
#ifdef NDEBUG
    sites = 7;
#else
    sites = 5;
#endif

    HubbardSystem s1(sites, true);
    DoubleOccupancy<HubbardSystem> DO1(s1);
    s1.H.construct();
    s1.H.diagonalizeUsingSymmetries();

    HubbardSystem s2(sites, true);
    DoubleOccupancy<HubbardSystem> DO2(s2);
    s2.H.construct();
    s2.H.diagonalizeUsingSymmetriesBySectors();

    startCPUTime = std::clock();
    s1.measure(1, DO1(0));
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for ensembleAverage of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s2.measureBS(1, DO2(0));
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for ensembleAverageBySectors of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;
}
