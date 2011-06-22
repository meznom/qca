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
    HubbardHamiltonian (const System& s_) 
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
                //t=1
                H += s.creator(2*i+spin) * s.annihilator(2* ((i+1)%s.N_sites) +spin);
                H += s.creator(2* ((i+1)%s.N_sites) +spin) * s.annihilator(2*i+spin);
            }
            H += 10 * s.creator(2*i) * s.annihilator(2*i) * //U=10 
                     s.creator(2*i+1) * s.annihilator(2*i+1);
        }
    }

    SMatrix& H;
    const System& s;
};

typedef BasicSystem BaseSystem;
class HubbardSystem : public BaseSystem
{
public:
    HubbardSystem (size_t N_sites_, bool exploitSymmetries_ = false) 
        : BaseSystem(2*N_sites_), N_sites(N_sites_), H(*this), 
          exploitSymmetries(exploitSymmetries_), ensembleAverage(*this)
    {}

    void construct ()
    {
        if (exploitSymmetries)
        {
            //TODO: technically there's a problem here, because N and S get
            //destroyed before basis gets destroyed
            basis.addSymmetryOperator(&N);
            basis.addSymmetryOperator(&S);
        }
        BaseSystem::construct();
        H.construct();
    }

    size_t N_sites;
    HubbardHamiltonian<HubbardSystem> H;
    bool exploitSymmetries;
    ParticleNumberSymmetryOperator N;
    SpinSymmetryOperator S;
    EnsembleAverage<HubbardSystem> ensembleAverage;
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
    s.construct();
    BOOST_CHECK (s.basis.getRanges().size() == 1);
    
    s.H.diagonalize();
    DVector eigenvalues1 = s.H.eigenvalues;
    SMatrix eigenvectors1 = s.H.eigenvectors;
    
    s.H.diagonalizeBlockWise();
    DVector eigenvalues2 = s.H.eigenvalues;
    SMatrix eigenvectors2 = s.H.eigenvectors;

    BOOST_CHECK (eigenvalues1 == eigenvalues2);
    BOOST_CHECK (eigenvectors1 == eigenvectors2);
}

BOOST_AUTO_TEST_CASE ( construct_system_with_symmetries )
{
    HubbardSystem s(4, true);
    s.construct();
    BOOST_CHECK (s.basis.getRanges().size() > 1);

    s.H.diagonalize();
    DVector eigenvalues1 = s.H.eigenvalues;
    DMatrix eigenvectors1 = EigenHelpers::sparseToDenseBlock(s.H.eigenvectors, 
                                                             0, 0, 
                                                             s.H.eigenvectors.cols(), 
                                                             s.H.eigenvectors.cols());
    std::vector<double> ev1;
    for (int i=0; i<eigenvalues1.size(); i++)
        ev1.push_back(eigenvalues1(i));
    
    s.H.diagonalizeBlockWise();
    DVector eigenvalues2 = s.H.eigenvalues;
    DMatrix eigenvectors2 = EigenHelpers::sparseToDenseBlock(s.H.eigenvectors, 
                                                             0, 0, 
                                                             s.H.eigenvectors.cols(), 
                                                             s.H.eigenvectors.cols());
    std::vector<double> ev2;
    for (int i=0; i<eigenvalues2.size(); i++)
        ev2.push_back(eigenvalues2(i));

    BOOST_CHECK (ev1.size() == ev2.size());

    //We only check if both methods yield the same eigenvalues. Checking
    //whether the eigenvectors are equivalent is too complicated.
    for (size_t i=0; i<ev1.size(); i++)
        BOOST_CHECK (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i])) != ev2.end());
}

BOOST_AUTO_TEST_CASE ( measure_double_occupancy )
{
    HubbardSystem s(4, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.construct();
    s.H.diagonalize();
    double doLowT1 = s.ensembleAverage(1000, DO(0));
    double doMidT1 = s.ensembleAverage(1, DO(0));
    double doHighT1 = s.ensembleAverage(0.00001, DO(0));

    s.H.diagonalizeBlockWise();
    double doLowT2 = s.ensembleAverage(1000, DO(0));
    double doMidT2 = s.ensembleAverage(1, DO(0));
    double doHighT2 = s.ensembleAverage(0.00001, DO(0));

    BOOST_CHECK (epsilonEqual(doLowT1, 0, 10E-2));
    BOOST_CHECK (epsilonEqual(doHighT1, 0.25, 10E-2));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT2));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT2));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT2));
}

BOOST_AUTO_TEST_CASE ( compare_performance_diagonalize_and_diagonalizeBlockWise )
{
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;
    size_t sites;
#ifdef NDEBUG
    sites = 5;
#else
    sites = 4;
#endif
    HubbardSystem s(sites, true);
    DoubleOccupancy<HubbardSystem> DO(s);


    startCPUTime = std::clock();
    s.construct();
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
    s.H.diagonalizeBlockWise();
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time for diagonalizeBlockWise of system with " << sites 
              << " sites: " << cpuTime << "s" << std::endl;
}
