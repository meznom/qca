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
        s.basis.applyMask(H);
    }

    SMatrix& H;
    const System& s;
};

typedef BasicSystem<Filter::SelectAll, Sorter::DontSort> BaseSystem;
class HubbardSystem : public BaseSystem
{
public:
    HubbardSystem (size_t N_sites_, bool exploitSymmetries_ = false) 
        : BaseSystem(2*N_sites_, Filter::SelectAll(), Sorter::DontSort()), 
          N_sites(N_sites_), H(*this), exploitSymmetries(exploitSymmetries_),
          ensembleAverage(*this)
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
    typedef std::pair<double, DVector> ValueType;
    EpsilonEqualPred (ValueType v_, double epsilon_ = 10E-10)
        : v(v_), epsilon(epsilon_)
    {}

    bool operator() (const ValueType& w) const
    {
        //if (fabs(v.first - w.first) > epsilon)
        //    return false;
        assert(v.second.size() == w.second.size());
        for (int i=0; i<v.second.size(); i++)
            //if (fabs(v.first * v.second(i) - w.first * w.second(i)) > epsilon)
            if (fabs(v.second(i) - w.second(i)) > epsilon)
                return false;
        return true;
    }

private:
    ValueType v;
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
    
    s.H.diagonalise();
    DVector eigenvalues1 = s.H.eigenvalues;
    SMatrix eigenvectors1 = s.H.eigenvectors;
    
    s.H.diagonaliseBlockWise();
    DVector eigenvalues2 = s.H.eigenvalues;
    SMatrix eigenvectors2 = s.H.eigenvectors;

    BOOST_CHECK (eigenvalues1 == eigenvalues2);
    BOOST_CHECK (eigenvectors1 == eigenvectors2);
}

BOOST_AUTO_TEST_CASE ( construct_system_with_symmetries )
{
    HubbardSystem s(3, true);
    s.construct();
    BOOST_CHECK (s.basis.getRanges().size() == 9);

    s.H.diagonalise();
    DVector eigenvalues1 = s.H.eigenvalues;
    DMatrix eigenvectors1 = EigenHelpers::sparseToDenseBlock(s.H.eigenvectors, 
                                                             0, 0, 
                                                             s.H.eigenvectors.cols(), 
                                                             s.H.eigenvectors.cols());
    std::vector< std::pair<double, DVector> > ev1;
    for (int i=0; i<eigenvalues1.size(); i++)
        ev1.push_back(std::pair<double,DVector>(eigenvalues1(i), eigenvectors1.col(i)));
    
    s.H.diagonaliseBlockWise();
    DVector eigenvalues2 = s.H.eigenvalues;
    DMatrix eigenvectors2 = EigenHelpers::sparseToDenseBlock(s.H.eigenvectors, 
                                                             0, 0, 
                                                             s.H.eigenvectors.cols(), 
                                                             s.H.eigenvectors.cols());
    std::vector< std::pair<double, DVector> > ev2;
    for (int i=0; i<eigenvalues2.size(); i++)
        ev2.push_back(std::pair<double,DVector>(eigenvalues2(i), eigenvectors2.col(i)));

    //this is only a rough check, it does not properly handle degeneracies
    for (size_t i=0; i<ev1.size(); i++)
    {
        if (epsilonEqual(ev1[i].first, -2))
            std::cerr << "-> " << i << "  " << ev1[i].first << std::endl << ev1[i].second << std::endl << std::endl;
        //BOOST_CHECK (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i], 10E-3)) != ev2.end());
        //if (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i], 10E-3)) == ev2.end())
        //{
        //    std::cerr << "-> " << i << std::endl << ev1[i].first << std::endl << std::endl << ev1[i].second << std::endl << std::endl;
        //    std::cerr << "--> " << s.H.H * ev1[i].second << std::endl << std::endl << std::endl;
        //    break;
        //}
    }

    std::cerr << std::endl << std::endl << "Other way around: " << std::endl << std::endl;
    for (size_t i=0; i<ev2.size(); i++)
    {
        if (epsilonEqual(ev2[i].first, -2))
            std::cerr << "-> " << i << "  " << ev2[i].first << std::endl << ev2[i].second << std::endl << std::endl;
        //BOOST_CHECK (std::find_if(ev2.begin(), ev2.end(), EpsilonEqualPred(ev1[i], 10E-3)) != ev2.end());
        //if (std::find_if(ev1.begin(), ev1.end(), EpsilonEqualPred(ev2[i], 10E-3)) == ev1.end())
        //{
        //    std::cerr << "-> " << i << std::endl << ev2[i].first << std::endl << std::endl << ev2[i].second << std::endl << std::endl;
        //    std::cerr << "--> " << s.H.H * ev2[i].second << std::endl << std::endl << std::endl;
        //    break;
        //}
    }
}

BOOST_AUTO_TEST_CASE ( measure_double_occupancy )
{
    HubbardSystem s(4, true);
    DoubleOccupancy<HubbardSystem> DO(s);
    s.construct();
    s.H.diagonalise();
    double doLowT1 = s.ensembleAverage(1000, DO(0));
    double doMidT1 = s.ensembleAverage(1, DO(0));
    double doHighT1 = s.ensembleAverage(0.00001, DO(0));

    s.H.diagonaliseBlockWise();
    double doLowT2 = s.ensembleAverage(1000, DO(0));
    double doMidT2 = s.ensembleAverage(1, DO(0));
    double doHighT2 = s.ensembleAverage(0.00001, DO(0));

    BOOST_CHECK (epsilonEqual(doLowT1, 0));
    BOOST_CHECK (epsilonEqual(doHighT1, 0.25, 10E-4));

    BOOST_CHECK (epsilonEqual(doLowT1, doLowT2));
    BOOST_CHECK (epsilonEqual(doMidT1, doMidT2));
    BOOST_CHECK (epsilonEqual(doHighT1, doHighT2));
}

//BOOST_AUTO_TEST_CASE ( compare_performance_diagonalise_and_diagonaliseBlockWise )
//{
//    std::clock_t startCPUTime, endCPUTime;
//    double cpuTime = 0;
//    HubbardSystem s(5, true);
//    DoubleOccupancy<HubbardSystem> DO(s);
//
//
//    startCPUTime = std::clock();
//    s.construct();
//    endCPUTime = std::clock();
//    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
//    std::cerr << "Time for construction: " << cpuTime << "s" << std::endl;
//
//    startCPUTime = std::clock();
//    s.H.diagonalise();
//    endCPUTime = std::clock();
//    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
//    std::cerr << "Time for diagonalise: " << cpuTime << "s" << std::endl;
//
//    startCPUTime = std::clock();
//    s.H.diagonaliseBlockWise();
//    endCPUTime = std::clock();
//    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
//    std::cerr << "Time for diagonaliseBlockWise: " << cpuTime << "s" << std::endl;
//}
