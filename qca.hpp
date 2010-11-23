#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"

class Hopping
{
public:
    Hopping (double t_, double td_)
    : t(t_), td(td_)
    {}

    double operator() (size_t i, size_t j) const
    {
        if (i == j)
            return 0;
        if (std::abs( static_cast<int>(i) - static_cast<int>(j) ) == 2)
            return td;
        return t;
    }

private:
    const double t, td;
};

class Coulomb
{
public:
    Coulomb (double V0_, double a_, double b_)
    : V0(V0_), a(a_), b(b_)
    {}

    double operator() (size_t i, size_t j) const
    {
        if (i == j)
            return V0;
        return 1 / distance(i,j);

    }

    double distance (size_t i_, size_t j_) const
    {
        const int i = static_cast<int>(i_);
        const int j = static_cast<int>(j_);
        const double deltaY = a * ( (i/2-j/2)%2 );
        assert( std::abs(deltaY-a) < 10E-20 || deltaY < 10E-20 );

        /*
         * 0 1
         * 3 2
         */
        //TODO: check if this is correct, can we optimise it?
        const double deltaX = 
            (a+b) * ( (i - j) / 4 ) + 
            a * ( ( ((i%4)%3==0)?0:1 ) - ( ((j%4)%3==0)?0:1 ) );

        //std::cerr << i << "   " << j << "    " << deltaX << "   " << deltaY << std::endl;

        return std::sqrt(deltaX*deltaX + deltaY*deltaY);
    }

private:
    const double V0, a, b;
};

class External
{
public:
    External (double Vext_)
    : Vext(Vext_)
    {}

    double operator() (size_t i) 
    {
        if (i==0)
            return Vext;
        if (i==3)
            return -Vext;
        return 0;
    }
private:
    double Vext;
};

namespace Filter {
class NElectronsPerPlaquet
{
public:
    NElectronsPerPlaquet (size_t N_, size_t plaquetSize_)
    : N(N_), plaquetSize(plaquetSize_)
    {}

    bool operator() (const State& s) const
    {
        if (s.count() != N*s.size()/plaquetSize)
            return false;
        for (size_t i=0; i<s.size(); i+=plaquetSize)
            if (s.count(i, i+plaquetSize) != N)
                return false;
        return true;
    }

private:
    size_t N, plaquetSize;
};
}; /* namespace Filter */

namespace Sorter {
class Bond
{
public:
    Bond ()
    : twoElectrons(2,4)
    {}

    bool operator() (const State& s1, const State& s2) const
    {
        /*
         * Sort by particle number first
         */
        if (s1.count() != s2.count())
            return s1.count() < s2.count();

        /*
         * Special sorting rule for half-filling, i.e. 2 electrons per plaquet.
         */
        //if (s1.count() != s1.size()/2)
        //    return false;
        if (twoElectrons(s1) && twoElectrons(s2))
        {
            for (size_t i=0; i<s1.size(); i+=4)
                if (stateNumber(s1,i) != stateNumber(s2,i))
                    return stateNumber(s1,i) < stateNumber(s2,i);
            return false;
        }
        /*
         * Enforce a definite ordering of twoElectron vs non-twoElectron per
         * plaquet systems
         */
        else
            return twoElectrons(s1) && !twoElectrons(s2);
    }

    int stateNumber (const State& s, int offset) const
    {
        const int o = offset;
        if (s[0+o] && s[1+o]) return 1;
        if (s[1+o] && s[2+o]) return 2;
        if (s[2+o] && s[3+o]) return 3;
        if (s[3+o] && s[0+o]) return 4;
        if (s[0+o] && s[2+o]) return 5;
        if (s[1+o] && s[3+o]) return 6;
        return 0;
    }

private:
    Filter::NElectronsPerPlaquet twoElectrons;
};
}; /* namespace Sorter */

template<class System>
class Polarisation
{
public:
    Polarisation (const System& s_)
    : s(s_)
    {}

    SMatrix operator() (size_t p) const
    {
        const size_t o = 4*p;
        return s.basis.applyMask(
            1.0/4.0 * ( s.n(o+1)+s.n(o+3) - s.n(o+0)-s.n(o+2) )
        );
    }

private:
    const System& s;
};

template<class System>
class CreatorAnnihilator
{
public:
    CreatorAnnihilator (const System& s_, size_t plaquetSize_)
    : s(s_), plaquetSize(plaquetSize_), cas(s.N_p * plaquetSize * plaquetSize), 
      zeroMatrix(s.basis.size(), s.basis.size())
    {
        //TODO: optimise - c_i a_j = (c_j a_i)^{\dag}
        for (size_t p=0; p<s.N_p; p++)
            for (size_t i=0; i<plaquetSize; i++)
                for (size_t j=0; j<plaquetSize; j++)
                    constructMatrix(plaquetSize*p + i, plaquetSize*p + j);
    }

    const SMatrix& operator() (size_t i, size_t j) const
    {
        // no interplaquet hopping => return a 0-matrix
        if (i/plaquetSize != j/plaquetSize)
            return zeroMatrix;
        return cas[I(i,j)];
    }


private:
    void constructMatrix (size_t i, size_t j)
    {
        SMatrix& m = cas[I(i,j)];
        m = SMatrix(s.basis.size(), s.basis.size());
        // we expect one entry per column
        m.reserve(s.basis.size());
        for (size_t col=0; col<s.basis.size(); col++)
        {
            m.startVec(col);
            if (i==j && s.basis(col)[i] == 1)
            {
                m.insertBack(col, col) = 1;
                continue;
            }
            if (s.basis(col)[i] == 1 || s.basis(col)[j] == 0)
                continue;
            State state(s.basis(col));
            state[i] = 1;
            state[j] = 0;
            const size_t row = s.basis(state);
            size_t sum = state.count(i,j); //TODO: is this correct?
            if (i<j) sum -= 1; //works, because for i<j we always have sum>=1
            const double sign = (sum%2==0)?1:-1; // probably faster than using (-1)^sum
            m.insertBack(row, col) = sign;
        }
        m.finalize();
    }

    size_t I (size_t i, size_t j) const
    {
        const size_t p = i/plaquetSize;
        assert(p == j/plaquetSize);
        const size_t ii = i%plaquetSize;
        const size_t jj = j%plaquetSize;
        return plaquetSize * plaquetSize * p + plaquetSize * ii + jj;
    }

    const System& s;
    const size_t plaquetSize;
    std::vector<SMatrix> cas;
    SMatrix zeroMatrix;
};

template<class System>
class QCAHamiltonian : public Hamiltonian<System>
{
public:
    QCAHamiltonian (const System& s_)
    : Hamiltonian<System>(s_), 
      t(1), td(0), V0(1000), a(1.0), b(3*a), Vext(0), 
      H(Hamiltonian<System>::H), s(Hamiltonian<System>::s)
    {}

    void construct() 
    {
        Hopping hopping(t, td);
        Coulomb coulomb(V0, a, b);
        External external(Vext);

        H.setZero();
        for (size_t i=0; i<s.N_sites; i++)
        {
            H += coulomb(i,i) * s.n_updown(i);
            H += external(i) * s.n(i);
            for (size_t j=i+1; j<s.N_sites; j++)
            {
                H += - hopping(i,j) * s.ca(i,j) - hopping(j,i) * s.ca(j,i);
                H += coulomb(i,j) * s.n(i) * s.n(j);
            }
        }

        s.basis.applyMask(H);
    }

    double t, td, V0, a, b, Vext;

    SMatrix& H;
    const System& s;
};

typedef MinimalSystem<Filter::NElectronsPerPlaquet, Sorter::Bond> QCABondBase;
class QCABond : public QCABondBase
{
public:
    QCABond (size_t N_p_)
    : QCABondBase(4*N_p_, Filter::NElectronsPerPlaquet(2,4), Sorter::Bond()), N_p(N_p_), 
      N_sites(4*N_p), ca(*this,4), H(*this), ensembleAverage(*this), P(*this)
    {}

    /*
    SMatrix ca (size_t i, size_t j) const
    {
        return creator(i) * annihilator(j);
    }
    */

    SMatrix n (size_t i) const
    {
        return ca(i,i);
    }

    SMatrix n_updown (size_t i) const
    {
        // return 0
        return SMatrix(basis.size(), basis.size());
    }

    size_t N_p, N_sites;
    CreatorAnnihilator<QCABond> ca;
    QCAHamiltonian<QCABond> H;
    EnsembleAverage<QCABond> ensembleAverage;
    Polarisation<QCABond> P;
};

namespace Filter {
class Spin
{
public:
    Spin (int S_)
    : S(S_)
    {}


    bool operator() (const State& s) const
    {
        //TODO
         return spin(s) == S;
    }

    
    int spin (const State& s) const
    {
        // N_down = N - N_up with N the total particle number
        // spin = N_up - N_down = 2 N_up - N
        const int N = s.count();
        int N_up=0;
        for (size_t i=0; i<s.size(); i+=2)
            N_up += s[i];
        return 2*N_up - N;
    }

private:
    const int S;
};

template<class Filter1, class Filter2>
class And
{
public:
    And (Filter1 f1_, Filter2 f2_)
    : f1(f1_), f2(f2_)
    {}

    bool operator() (const State& s) const
    {
        return f1(s) && f2(s);
    }
private:
    Filter1 f1;
    Filter2 f2;
};
}; /* namespace Filter */

namespace Sorter {
class ParticleNumberAndSpin
{
public:
    bool operator() (const State& s1, const State& s2) const
    {
        for (size_t i=0; i<s1.size(); i+=8)
            if (s1.count(i,i+8) != s2.count(i,i+8))
                return s1.count(i,i+8) < s2.count(i,i+8);

        return spin(s1) < spin(s2);
    }

    int spin (const State& s) const
    {
        // N_down = N - N_up with N the total particle number
        // spin = N_up - N_down = 2 N_up - N
        const int N = s.count();
        int N_up=0;
        for (size_t i=0; i<s.size(); i+=2)
            N_up += s[i];
        return 2*N_up - N;
    }
};
}; /* namespace Sorter */

typedef MinimalSystem<Filter::NElectronsPerPlaquet, Sorter::ParticleNumberAndSpin> QCAQuarterFillingBase;
class QCAQuarterFilling : public QCAQuarterFillingBase
{
public:
    QCAQuarterFilling (size_t N_p_)
    : QCAQuarterFillingBase(8*N_p_, Filter::NElectronsPerPlaquet(2,8), Sorter::ParticleNumberAndSpin()), 
      N_p(N_p_), N_sites(4*N_p), creatorAnnihilator(*this,8), H(*this), 
      ensembleAverage(*this), P(*this)
    {}

    size_t I (size_t i, Spin s) const
    {
        return 2*i + s;
    }

    /*
    SMatrix ca (size_t i, Spin s_i, size_t j, Spin s_j) const
    {
        return creator(I(i, s_i))*annihilator(I(j, s_j));
    }
    */

    SMatrix ca (size_t i, Spin s_i, size_t j, Spin s_j) const
    {
        return creatorAnnihilator(I(i, s_i), I(j, s_j));
    }

    SMatrix ca (size_t i, size_t j) const
    {
        return ca(i,UP,j,UP) + ca(i,DOWN,j,DOWN);
    }

    SMatrix n (size_t i, Spin s) const
    {
        return ca(i,s,i,s);
    }

    SMatrix n (size_t i) const
    {
        return n(i,UP) + n(i,DOWN);
    }

    SMatrix n_updown (size_t i) const
    {
        return n(i,UP) * n(i,DOWN);
    }

    size_t N_p, N_sites;
    CreatorAnnihilator<QCAQuarterFilling> creatorAnnihilator;
    QCAHamiltonian<QCAQuarterFilling> H;
    EnsembleAverage<QCAQuarterFilling> ensembleAverage;
    Polarisation<QCAQuarterFilling> P;
};


#endif // __QCA_HPP__
