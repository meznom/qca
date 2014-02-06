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
    SMatrix& H;
    System& s;

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
};

class HubbardSystem
{
public:
    typedef HubbardSystem Self;

    size_t N_sites;
    size_t N_orbitals;
    
    Basis basis;
    Creator<Self> creator;
    Annihilator<Self> annihilator;
    HubbardHamiltonian<Self> H;
    EnsembleAverage<Self> measure;
    
    ParticleNumberSymmetryOperator N;
    SpinSymmetryOperator S;
    bool exploitSymmetries;
    double t, U;

public:
    HubbardSystem (size_t N_sites_, bool exploitSymmetries_ = false) 
        : N_sites(N_sites_), N_orbitals(2*N_sites_), 
          creator(*this), annihilator(*this), H(*this), measure(*this),
          exploitSymmetries(exploitSymmetries_), t(1), U(10)
    {
        if (exploitSymmetries)
        {
            basis.addSymmetryOperator(&N);
            basis.addSymmetryOperator(&S);
        }
        construct();
    }
    
    void construct ()
    {
        basis.construct(N_orbitals);
        creator.construct();
        annihilator.construct();
    }

    void update()
    {
        H.construct();
        H.diagonalize();
    }
};

template<class System>
class DoubleOccupancy
{
private:
    const System& s;

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
    
    s.H.diagonalize();
    DVector eigenvalues = s.H.eigenvalues();
    auto& eigenvectors = s.H.eigenvectorsBySector();

    BOOST_CHECK (eigenvalues.size() == 16);
    BOOST_CHECK (eigenvectors.size() == 1);
    BOOST_CHECK (eigenvectors[0].cols() == 16);
}

BOOST_AUTO_TEST_CASE ( construct_system_with_symmetries )
{
    HubbardSystem s(4, true);
    s.H.construct();
    BOOST_CHECK (s.basis.getRanges().size() > 1);

    s.H.diagonalize();
    DVector eigenvalues = s.H.eigenvalues();
    auto& eigenvectors = s.H.eigenvectorsBySector();
    
    BOOST_CHECK (eigenvalues.size() == 256);
    BOOST_CHECK (eigenvectors.size() == 25);
    int n = 0;
    for (auto i=eigenvectors.begin(); i!=eigenvectors.end(); i++)
        n += i->cols();
    BOOST_CHECK (n == 256);
}

BOOST_AUTO_TEST_CASE ( construct_and_diagonalize_system_multiple_times )
{
    HubbardSystem s(4, true);
    s.H.construct();
    size_t n_b1 = s.basis.size();
    size_t n_r1 = s.basis.getRanges().size();
    
    s.H.diagonalize();
    DVector ev1 = s.H.eigenvalues();

    s.H.construct();
    size_t n_b2 = s.basis.size();
    size_t n_r2 = s.basis.getRanges().size();
    BOOST_CHECK (n_b1 == n_b2);
    BOOST_CHECK (n_r1 == n_r2);
   
    s.H.diagonalize();
    DVector ev2 = s.H.eigenvalues();

    BOOST_CHECK (ev1 == ev2);
}

BOOST_AUTO_TEST_CASE ( test_by_sector_eigenvalue_and_eigenvectors_accessors )
{
    HubbardSystem s(4, true);
    s.H.construct();
   
    s.H.diagonalize();
    auto eigenvalues1 = s.H.eigenvaluesBySector();
    auto energies1 = s.H.eigenvalues();

    BOOST_CHECK (eigenvalues1.size() > 1);
    BOOST_CHECK (eigenvalues1[0].size() > 0);
    BOOST_CHECK (energies1.size() == 256);

    DVector ev(256);
    int k = 0;
    for (size_t j=0; j<eigenvalues1.size(); j++)
        for (int i=0; i<eigenvalues1[j].size(); i++)
        {
            ev(k) = eigenvalues1[j](i);
            k++;
        }
    BOOST_CHECK (k == 256);
    BOOST_CHECK (energies1 == ev);

    // test that eigenvectors get properly reset and reconstructed
    s.H.diagonalize();
    auto energies2 = s.H.eigenvalues();
    BOOST_CHECK (energies2.size() == 256);

    s.t = 1000;
    s.U = 20000;
    s.H.construct();
    s.H.diagonalize();
    auto energies3 = s.H.eigenvalues();

    double diff = 0;
    for (int i=0; i<256; i++)
        diff += std::abs(energies2(i) - energies3(i));
    BOOST_CHECK ( diff > 1 );
}

BOOST_AUTO_TEST_CASE ( measure_double_occupancy )
{
    HubbardSystem s(4, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.H.construct();
    
    s.H.diagonalize();
    double doLowT = s.measure(1000, DO(0));
    double doMidT = s.measure(1, DO(0));
    double doHighT = s.measure(0.00001, DO(0));

    BOOST_CHECK (epsilonEqual(doLowT, 0, 10E-2));
    BOOST_CHECK (epsilonEqual(doMidT, 0.00453281, 10E-6));
    BOOST_CHECK (epsilonEqual(doHighT, 0.25, 10E-2));
}

BOOST_AUTO_TEST_CASE ( performance_of_diagonalization )
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
    s.H.diagonalize();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalize of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;

    startCPUTime = std::clock();
    s.measure(1, DO(0));
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for ensembleAverage of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;
}
