#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"
#include "utilities.hpp"

class Hopping
{
public:
    Hopping (double t_, double td_, double ti_)
    : t(t_), td(td_), ti(ti_)
    {}

    double operator() (size_t i, size_t j) const
    {
        /*
         * Same plaquet
         */
        if (i/4 == j/4)
        {
            if (std::abs( static_cast<int>(i) - static_cast<int>(j) ) == 2)
                return td;
            else if (i != j)
                return t;
        }
        /*
         * Neighbouring plaquets
         */
        else if (std::abs( static_cast<int>(i/4) - static_cast<int>(j/4) ) == 1)
        {
            //TODO: this is untested
            const size_t l = std::min(i, j); //left plaquet
            const size_t r = std::max(i, j); //right plaquet
            if ( (l%4 == 1 && r%4 == 0) || (l%4 == 2 && r%4 == 3) )
                return ti;
        }
        return 0;
    }

private:
    const double t, td, ti;
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
         * 0 1   4 5   ...
         * 3 2   7 6   ...
         *
         * verified for 1, 2 and 3 plaquets, so seems to work correctly
         */
        const double deltaX = 
            (a+b) * (i/4 - j/4) + 
            a * ( ( ((i%4)%3==0)?0:1 ) - ( ((j%4)%3==0)?0:1 ) );

        //std::cerr << i << "   " << j << "    " << deltaX << "   " << deltaY << std::endl;

        return std::sqrt(deltaX*deltaX + deltaY*deltaY);
    }

private:
    const double V0, a, b;
};

template<class ParameterContainer>
class ExternalPlain
{
public:
    ExternalPlain (const ParameterContainer& c)
    : Vext(c.Vext)
    {}

    double operator() (size_t i) const
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

template<class ParameterContainer>
class ExternalDeadPlaquet
{
public:
    ExternalDeadPlaquet (const ParameterContainer& c)
    : coulomb(c.V0, c.a, c.b), P(c.Pext)
    {}

    double operator() (size_t i) const
    {
        /*
         * Physically the dead plaquet sits to the left of the linear chain
         * system, at -4,-3,-2,-1. To use our Coulomb distance method we shift
         * the whole system one plaquet to the right (and thus all sites are
         * positive).
         */
        double V = 0;
        for (int j=1; j<4; j+=2)
            V += (P+1)/2 * 1/coulomb.distance(j,i+4);
        for (int j=0; j<4; j+=2)
            V += (1-P)/2 * 1/coulomb.distance(j,i+4);
        return V;
    }
private:
    const Coulomb coulomb;
    const double P;
};

template<class System, template <typename> class External>
class QcaHamiltonian : public Hamiltonian<System>
{
public:
    QcaHamiltonian (const System& s_)
    : Hamiltonian<System>(s_), 
      t(1), td(0), ti(0), V0(1000), a(1.0), b(3*a), Vext(0), Pext(0), mu(0), 
      H(Hamiltonian<System>::H), s(Hamiltonian<System>::s)
    {}

    void construct() 
    {
        Hopping hopping(t, td, ti);
        Coulomb coulomb(V0, a, b);
        External<QcaHamiltonian> external(*this);

        H.setZero();
        for (size_t i=0; i<s.N_sites; i++)
        {
            H += coulomb(i,i) * s.n_updown(i);
            //H += (externalPlain(i) + externalDP(i) + mu) * s.n(i);
            H += (external(i) + mu) * s.n(i);
            for (size_t j=i+1; j<s.N_sites; j++)
            {
                H += - hopping(i,j) * s.ca(i,j) - hopping(j,i) * s.ca(j,i);
                H += coulomb(i,j) * s.n(i) * s.n(j);
            }
        }

        std::cerr << "Number of non-zero elements on H: " << H.nonZeros() << std::endl;

        s.basis.applyMask(H);
    }

    double t, td, ti, V0, a, b, Vext, Pext, mu;

    SMatrix& H;
    const System& s;
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
            1.0/2.0 * ( s.n(o+1)+s.n(o+3) - s.n(o+0)-s.n(o+2) )
        );
    }

private:
    const System& s;
};

template<class System>
class ParticleNumber
{
public:
    ParticleNumber (const System& s_)
    : s(s_)
    {}

    SMatrix operator() (size_t p) const
    {
        const size_t o = 4*p;
        return s.basis.applyMask(
            s.n(0) + s.n(1) + s.n(2) + s.n(3)
        );
    }

    SMatrix operator() () const
    {
        SMatrix N;
        for (size_t i=0; i<s.N_sites; i++)
            N += s.n(i);
        return N;
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

typedef MinimalSystem<Filter::NElectronsPerPlaquet, Sorter::Bond> QcaBondBase;
template<template <typename> class External>
class QcaBond : public QcaBondBase
{
public:
    QcaBond (size_t N_p_)
    : QcaBondBase(4*N_p_, Filter::NElectronsPerPlaquet(2,4), Sorter::Bond()), N_p(N_p_), 
      N_sites(4*N_p), ca(*this,4), H(*this), ensembleAverage(*this), P(*this), N(*this)
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
    CreatorAnnihilator<QcaBond> ca;
    QcaHamiltonian<QcaBond, External> H;
    EnsembleAverage<QcaBond> ensembleAverage;
    Polarisation<QcaBond> P;
    ParticleNumber<QcaBond> N;
};

typedef MinimalSystem<Filter::NElectronsPerPlaquet, Sorter::ParticleNumberAndSpin> QcaQuarterFillingBase;
template<template <typename> class External>
class QcaQuarterFilling : public QcaQuarterFillingBase
{
public:
    QcaQuarterFilling (size_t N_p_)
    : QcaQuarterFillingBase(8*N_p_, Filter::NElectronsPerPlaquet(2,8), Sorter::ParticleNumberAndSpin()), 
      N_p(N_p_), N_sites(4*N_p), creatorAnnihilator(*this,8), H(*this), 
      ensembleAverage(*this), P(*this), N(*this)
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
    CreatorAnnihilator<QcaQuarterFilling> creatorAnnihilator;
    QcaHamiltonian<QcaQuarterFilling, External> H;
    EnsembleAverage<QcaQuarterFilling> ensembleAverage;
    Polarisation<QcaQuarterFilling> P;
    ParticleNumber<QcaQuarterFilling> N;
};

typedef BasicSystem<Filter::SelectAll, Sorter::ParticleNumberAndSpin> QcaGrandCanonicalBase;
template<template <typename> class External>
class QcaGrandCanonical : public QcaGrandCanonicalBase
{
public:
    QcaGrandCanonical (size_t N_p_)
    : QcaGrandCanonicalBase(8*N_p_, Filter::SelectAll(), Sorter::ParticleNumberAndSpin()),
      N_p(N_p_), N_sites(4*N_p_), H(*this), ensembleAverage(*this), P(*this), N(*this)
    {}
    
    size_t I (size_t i, Spin s) const
    {
        return 2*i + s;
    }

    SMatrix ca (size_t i, Spin s_i, size_t j, Spin s_j) const
    {
        return creator(I(i, s_i))*annihilator(I(j, s_j));
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
    QcaHamiltonian<QcaGrandCanonical, External> H;
    EnsembleAverage<QcaGrandCanonical> ensembleAverage;
    Polarisation<QcaGrandCanonical> P;
    ParticleNumber<QcaGrandCanonical> N;
};

template<class QcaSystem>
class DQcaGeneric : public QcaSystem
{
public:
    DQcaGeneric (Description desc_)
    : QcaSystem (desc_["N_p"]), desc(desc_)
    {
        QcaSystem::H.t = desc["t"].get<double>(1.0);
        QcaSystem::H.td = desc["td"].get<double>(0); 
        QcaSystem::H.ti = desc["ti"].get<double>(0); 
        QcaSystem::H.a = desc["a"].get<double>(1.0); 
        QcaSystem::H.b = desc["b"].get<double>(3);
        QcaSystem::H.Vext = desc["Vext"].get<double>(0);
        QcaSystem::H.Pext = desc["Pext"].get<double>(0);
        QcaSystem::H.V0 = desc["V0"].get<double>(1000); 
        QcaSystem::H.mu = desc["mu"].get<double>(0);
    }

private:
    Description desc;
};

/*
 * Useful typedefs
 */
typedef QcaBond<ExternalPlain> QcaBondPlain;
typedef QcaBond<ExternalDeadPlaquet> QcaBondDeadPlaquet;
typedef QcaQuarterFilling<ExternalPlain> QcaQuarterFillingPlain;
typedef QcaQuarterFilling<ExternalDeadPlaquet> QcaQuarterFillingDeadPlaquet;
typedef QcaGrandCanonical<ExternalPlain> QcaGrandCanonicalPlain;
typedef QcaGrandCanonical<ExternalDeadPlaquet> QcaGrandCanonicalDeadPlaquet;

typedef DQcaGeneric<QcaBond<ExternalPlain> > DQcaBondPlain;
typedef DQcaGeneric<QcaBond<ExternalDeadPlaquet> > DQcaBondDeadPlaquet;
typedef DQcaGeneric<QcaQuarterFilling<ExternalPlain> > DQcaQuarterFillingPlain;
typedef DQcaGeneric<QcaQuarterFilling<ExternalDeadPlaquet> > DQcaQuarterFillingDeadPlaquet;
typedef DQcaGeneric<QcaGrandCanonical<ExternalPlain> > DQcaGrandCanonicalPlain;
typedef DQcaGeneric<QcaGrandCanonical<ExternalDeadPlaquet> > DQcaGrandCanonicalDeadPlaquet;

/*
class DQcaBond : public QcaBondPlain
{
public:
    DQcaBond (Description desc_)
    : QcaBondPlain(desc_["N_p"]), desc(desc_)
    {
        H.t = desc["t"].get<double>(1.0);
        H.td = desc["td"].get<double>(0); 
        H.ti = 0;
        H.a = desc["a"].get<double>(1.0); 
        H.b = desc["b"].get<double>(3);
        H.Vext = desc["Vext"].get<double>(0);
        H.Pext = desc["Pext"].get<double>(0);
        H.V0 = 0; 
        H.mu = 0;
    }

private:
    Description desc;
};

class DQcaQuarterFilling: public QcaQuarterFilling
{
public:
    DQcaQuarterFilling (Description desc_)
    : QcaQuarterFilling(desc_["N_p"]), desc(desc_)
    {
        H.t = desc["t"].get<double>(1.0);
        H.td = desc["td"].get<double>(0); 
        H.ti = 0;
        H.a = desc["a"].get<double>(1.0); 
        H.b = desc["b"].get<double>(3);
        H.Vext = desc["Vext"].get<double>(0);
        H.Pext = desc["Pext"].get<double>(0);
        H.V0 = desc["V0"].get<double>(1000); 
        H.mu  = 0;
    }

private:
    Description desc;
};

class DQcaGrandCanonical: public QcaGrandCanonical
{
public:
    DQcaGrandCanonical (Description desc_)
    : QcaGrandCanonical (desc_["N_p"]), desc(desc_)
    {
        H.t = desc["t"].get<double>(1.0);
        H.td = desc["td"].get<double>(0); 
        H.ti = desc["ti"].get<double>(0); 
        H.a = desc["a"].get<double>(1.0); 
        H.b = desc["b"].get<double>(3);
        H.Vext = desc["Vext"].get<double>(0);
        H.Pext = desc["Pext"].get<double>(0);
        H.V0 = desc["V0"].get<double>(1000); 
        H.mu = desc["mu"].get<double>(0);
    }

private:
    Description desc;
};
*/

#endif // __QCA_HPP__
