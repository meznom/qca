#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE system test
#include <boost/test/unit_test.hpp>

#include "basis.hpp"
#include "system.hpp"
#include "qca.hpp"

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
        for (int spin=0; spin<2; spin++)
            for (size_t i=0; i<s.N_sites; i++)
            {
                H += s.creator(2*i+spin) * s.annihilator(2* ((i+1)%s.N_sites) +spin);
                H += s.creator(2* ((i+1)%s.N_sites) +spin) * s.annihilator(2*i+spin);
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
    HubbardSystem (size_t N_sites_) 
        : BaseSystem(2*N_sites_, Filter::SelectAll(), Sorter::DontSort()), 
          N_sites(N_sites_), H(*this)
    {}

    void construct ()
    {
        BaseSystem::construct();
        H.construct();
    }

    size_t N_sites;
    HubbardHamiltonian<HubbardSystem> H;
};


BOOST_AUTO_TEST_CASE ( construct_system )
{
    //QcaGrandCanonicalDeadPlaquet s(1);
    HubbardSystem s(2);
    s.construct();
    std::cerr << "number of ranges: " << s.basis.getRanges().size() << std::endl;
    std::cerr << s.basis.getRanges().front().a << "  " << s.basis.getRanges().front().b << std::endl;
    s.H.diagonalise();
    std::cerr << "diagonalise: " << s.H.eigenvalues << std::endl;
    std::cerr << "diagonalise: " << s.H.eigenvectors << std::endl;
    s.H.diagonaliseBlockWise();
    std::cerr << "diagonaliseBlockWise: " << s.H.eigenvalues << std::endl;
    std::cerr << "diagonaliseBlockWise: " << s.H.eigenvectors << std::endl;
}
