#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"
#include "utilities.hpp"
#include <limits>

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

        H = SMatrix(s.basis.size(), s.basis.size());
        H.setZero();
        for (size_t i=0; i<s.N_sites; i++)
        {
            H += coulomb(i,i) * s.n_updown(i);
            H += (external(i) + mu) * s.n(i);
            for (size_t j=i+1; j<s.N_sites; j++)
            {
                H += - hopping(i,j) * s.ca(i,j) - hopping(j,i) * s.ca(j,i);
                H += coulomb(i,j) * s.n(i) * s.n(j);
            }
        }
    }

    double t, td, ti, V0, a, b, Vext, Pext, mu;

    SMatrix& H;
    const System& s;
};

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
        return 1.0/2.0 * ( s.n(o+1)+s.n(o+3) - s.n(o+0)-s.n(o+2) );
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
        return s.n(o+0) + s.n(o+1) + s.n(o+2) + s.n(o+3);
    }

    SMatrix operator() () const
    {
        SMatrix N(s.basis.size(), s.basis.size());
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
    : s(s_), plaquetSize(plaquetSize_), cas(s.N_p * plaquetSize * plaquetSize) 
    {}

    void construct ()
    {
        zeroMatrix = SMatrix(s.basis.size(), s.basis.size());
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

class ParticleNumberPerPlaquetSymmetryOperator : public SymmetryOperator
{
public:
    ParticleNumberPerPlaquetSymmetryOperator (size_t plaquetSize_ = 8)
        : plaquetSize(plaquetSize_)
    {}

    virtual int operator() (const State& s) const
    {
        assert(s.size()/plaquetSize <= static_cast<size_t>(std::numeric_limits<int>::digits10));
        int N = 0;
        int multiplier = 1;
        for (size_t i=0; i<s.size(); i+=plaquetSize)
        {
            N += multiplier * s.count(i, i+plaquetSize);
            multiplier *= 10;
        }
        return N;
    }

    int valueForNElectronsPerPlaquet (int N, size_t N_p) const
    {
        assert(N_p <= static_cast<size_t>(std::numeric_limits<int>::digits10));
        int value = 0;
        int multiplier = 1;
        for (size_t i=0; i<N_p; i++)
        {
            value += N*multiplier;
            multiplier *= 10;
        }
        return value;
    }

private:
    size_t plaquetSize;
};

template<template <typename> class External>
class QcaBond
{
public:
    QcaBond (size_t N_p_)
        : basis(plaquetSize*N_p_), N_p(N_p_), N_sites(plaquetSize*N_p), 
          ca(*this, plaquetSize), H(*this), ensembleAverage(*this), P(*this), 
          N(*this), PPSO(plaquetSize)
    {}

    void construct ()
    {
        basis.addSymmetryOperator(&PPSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(2,N_p);
        basis.setFilter(constructSector(filterValue));
        basis.construct();
        ca.construct();
        H.construct();
    }

    SMatrix n (size_t i) const
    {
        return ca(i,i);
    }

    SMatrix n_updown (size_t i) const
    {
        // return 0
        return SMatrix(basis.size(), basis.size());
    }

    enum {plaquetSize=4};
    Basis basis;
    size_t N_p, N_sites;
    CreatorAnnihilator<QcaBond> ca;
    QcaHamiltonian<QcaBond, External> H;
    EnsembleAverage<QcaBond> ensembleAverage;
    Polarisation<QcaBond> P;
    ParticleNumber<QcaBond> N;
    ParticleNumberPerPlaquetSymmetryOperator PPSO;

private:
};

template<template <typename> class External>
class QcaQuarterFilling
{
public:
    QcaQuarterFilling (size_t N_p_)
        : basis(plaquetSize*N_p_), N_p(N_p_), N_sites(4*N_p), 
          creatorAnnihilator(*this, plaquetSize), H(*this), 
          ensembleAverage(*this), P(*this), N(*this), PPSO(plaquetSize)
    {}

    void construct ()
    {
        basis.addSymmetryOperator(&PPSO);
        basis.addSymmetryOperator(&SSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(2,N_p);
        basis.setFilter(constructSector(filterValue));
        basis.construct();
        creatorAnnihilator.construct();
        H.construct();
    }

    size_t I (size_t i, Spin s) const
    {
        return 2*i + s;
    }

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

    enum {plaquetSize=8};
    Basis basis;
    size_t N_p, N_sites;
    CreatorAnnihilator<QcaQuarterFilling> creatorAnnihilator;
    QcaHamiltonian<QcaQuarterFilling, External> H;
    EnsembleAverage<QcaQuarterFilling> ensembleAverage;
    Polarisation<QcaQuarterFilling> P;
    ParticleNumber<QcaQuarterFilling> N;
    ParticleNumberPerPlaquetSymmetryOperator PPSO;
    SpinSymmetryOperator SSO;
};

template<template <typename> class External>
class QcaGrandCanonical
{
public:
    QcaGrandCanonical (size_t N_p_)
        : basis(plaquetSize*N_p_), N_p(N_p_), N_sites(4*N_p_), creator(*this),
          annihilator(*this), H(*this), ensembleAverage(*this), P(*this), 
          N(*this)
    {}

    void construct ()
    {
        basis.addSymmetryOperator(&PSO);
        basis.addSymmetryOperator(&SSO);
        basis.construct();
        creator.construct();
        annihilator.construct();
        H.construct();
    }
    
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

    enum {plaquetSize=8};
    Basis basis;
    size_t N_p, N_sites;
    Creator<QcaGrandCanonical> creator;
    Annihilator<QcaGrandCanonical> annihilator;
    QcaHamiltonian<QcaGrandCanonical, External> H;
    EnsembleAverage<QcaGrandCanonical> ensembleAverage;
    Polarisation<QcaGrandCanonical> P;
    ParticleNumber<QcaGrandCanonical> N;
    ParticleNumberSymmetryOperator PSO;
    SpinSymmetryOperator SSO;
};

template<class QcaSystem>
class DQcaGeneric : public QcaSystem
{
public:
    DQcaGeneric (OptionSection os)
    : QcaSystem (os["p"])
    {
        setParameters(os);
    }

    void setParameters (OptionSection os)
    {
        QcaSystem::H.t = os["t"].get<double>(1.0);
        QcaSystem::H.td = os["td"].get<double>(0); 
        QcaSystem::H.ti = os["ti"].get<double>(0); 
        QcaSystem::H.a = os["a"].get<double>(1.0); 
        QcaSystem::H.b = os["b"].get<double>(3);
        QcaSystem::H.Vext = os["Vext"].get<double>(0);
        QcaSystem::H.Pext = os["Pext"].get<double>(0);
        QcaSystem::H.V0 = os["V0"].get<double>(1000); 
        QcaSystem::H.mu = os["mu"].get<double>(0);
    }
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

#endif // __QCA_HPP__
