#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"
#include "utilities.hpp"
#include <limits>

const double QCA_ELEMENTARY_CHARGE = 1.602176565E-19;
const double QCA_EPSILON_0 = 8.8541878176E-12;
const double QCA_NATURAL_EPSILON_R = QCA_ELEMENTARY_CHARGE / (4*M_PI*QCA_EPSILON_0*1e-9);

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
    Coulomb (double V0_, double a_, double b_, double lambdaD_ = 0, 
             double epsilonr_ = QCA_NATURAL_EPSILON_R, 
             double epsilon0_ = QCA_EPSILON_0)
    : V0(V0_), a(a_), b(b_), epsilon0(epsilon0_), epsilonr(epsilonr_), lambdaD(lambdaD_)
    {}

    double operator() (size_t i, size_t j) const
    {
        if (i == j)
            return V0;
        const double r = distance(i,j);
        if (lambdaD == 0)
            return QCA_ELEMENTARY_CHARGE / 
                   (4*M_PI * epsilon0 * epsilonr * r * 1e-9);
        else
            return QCA_ELEMENTARY_CHARGE * exp(- r / lambdaD) / 
                   (4*M_PI * epsilon0 * epsilonr * r * 1e-9);
    }

    double distance (size_t i_, size_t j_) const
    {
        const int i = static_cast<int>(i_);
        const int j = static_cast<int>(j_);
        
        /*
         * We compute deltaY and deltaX in units of a. This increases accuracy
         * for the case that a becomes very small. a and b are expected to be
         * roughly of the same order. Thus b/a~1 and deltaX~1 and deltaY~1. In
         * contrast we might have a~10^-10.
         */
        const double deltaY = (i/2-j/2)%2;
        assert( std::abs(deltaY-1) < 10E-20 || deltaY < 10E-20 );

        /*
         * 0 1   4 5   ...
         * 3 2   7 6   ...
         *
         * verified for 1, 2 and 3 plaquets, so seems to work correctly
         */
        const double deltaX = 
            (1+b/a) * (i/4 - j/4) + 
            ( ((i%4)%3==0)?0:1 ) - ( ((j%4)%3==0)?0:1 );

        return a * std::sqrt(deltaX*deltaX + deltaY*deltaY);
    }

private:
    const double V0, a, b, epsilon0, epsilonr, lambdaD;
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
    : coulomb(c.V0, c.a, c.b, c.lambdaD, c.epsilonr, c.epsilon0), 
      P(c.Pext), electronsPerPlaquet(c.electronsPerPlaquet), q(c.q)
    {}

    double operator() (size_t i) const
    {
        /*
         * Physically the dead plaquet sits to the left of the linear chain
         * system, at -4,-3,-2,-1. To use our Coulomb distance method we shift
         * the whole system one plaquet to the right (and thus all sites are
         * positive).
         */
        assert (electronsPerPlaquet == 2 || electronsPerPlaquet == 6);
        
        // Coulomb term: V_ij (n_i - q) (n_j - q)
        // here we are calculating V_ij (n_i - q) where i runs over all sites of
        // the dead plaquet
        double n = 0; // electrons per site
        if (electronsPerPlaquet == 6) n = 1;
        n -= q; //compensation charge
        double V = 0;
        for (int j=0; j<4; j++)
            V += n * coulomb(j,i+4);
        // put on two extra electrons according to the set polarization
        for (int j=1; j<4; j+=2)
            V += (P+1)/2 * coulomb(j,i+4);
        for (int j=0; j<4; j+=2)
            V += (1-P)/2 * coulomb(j,i+4);
        return V;
    }
private:
    const Coulomb coulomb;
    const double P;
    const size_t electronsPerPlaquet;
    const double q;
};

template<class System>
class QcaHamiltonian : public Hamiltonian<System>
{
public:
    QcaHamiltonian (const System& s_)
    : Hamiltonian<System>(s_), H(Hamiltonian<System>::H), 
      s(Hamiltonian<System>::s)
    {}

    void construct() 
    {
        Hopping hopping(s.t, s.td, s.ti);
        Coulomb coulomb(s.V0, s.a, s.b, s.lambdaD, s.epsilonr, s.epsilon0);
        typename System::External external(s);

        if (I.cols() == 0) constructIdentityMatrix();
        H = SMatrix(s.basis.size(), s.basis.size());
        H.setZero();
        for (size_t i=0; i<s.N_sites; i++)
        {
            /*
             * Eigen seems to truncate very small values in Sparse matrices. The
             * epsilon value used for the truncation is defined in Eigen's
             * NumTraits. To be safe we check that our values are bigger than
             * this threshold. 
             */
            assert(coulomb(i,i)==0 || 
                   std::fabs(coulomb(i,i))>NumTraits<double>::dummy_precision());
            assert((external(i)+s.mu)==0 || 
                   std::fabs((external(i)+s.mu))> NumTraits<double>::dummy_precision());
            
            H += coulomb(i,i) * s.n_updown(i);
            H += (external(i) - s.mu) * s.n(i);
            for (size_t j=i+1; j<s.N_sites; j++)
            {
                assert(hopping(i,j)==0 || 
                       std::fabs(hopping(i,j))>NumTraits<double>::dummy_precision());
                assert(hopping(j,i)==0 || 
                       std::fabs(hopping(j,i))>NumTraits<double>::dummy_precision());
                assert(coulomb(i,j)==0 || 
                       std::fabs(coulomb(i,j))>NumTraits<double>::dummy_precision());
                
                H += - hopping(i,j) * s.ca(i,j) - hopping(j,i) * s.ca(j,i);
                // V_ij (n_i - q) (n_j - q)
                H += coulomb(i,j) * (s.n(i) * s.n(j) - s.q * s.n(i) - 
                                     s.q * s.n(j) - s.q * s.q * I);
            }
        }
    }

private:
    void constructIdentityMatrix ()
    {
        I = SMatrix(s.basis.size(), s.basis.size());
        I.setZero();
        for (int i=0; i<I.cols(); i++)
        {
            I.startVec(i);
            I.insertBack(i,i) = 1;
        }
        I.finalize();
    }

    SMatrix& H;
    SMatrix I;
    const System& s;
};

template<class System>
class Polarization
{
public:
    Polarization (const System& s_)
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

template<class QcaSystem>
class QcaCommon
{
private:
    typedef QcaSystem S;
    S& s;

public:
    QcaCommon (QcaSystem& s_, size_t N_p_, size_t electronsPerPlaquet_ = 2)
        : s(s_), N_p(N_p_), N_sites(4*N_p), electronsPerPlaquet(electronsPerPlaquet_), 
          H(s), ensembleAverage(s), P(s), N(s), 
          t(1), td(0), ti(0), V0(1000), a(1.0), b(3*a), Vext(0), Pext(0), mu(0),
          epsilonr(QCA_NATURAL_EPSILON_R), lambdaD(0), 
          epsilon0(QCA_EPSILON_0), q(0)
    {
        assert(electronsPerPlaquet == 2 || electronsPerPlaquet == 6);
    }

    void update ()
    {
        H.construct();
        H.diagonalizeUsingSymmetriesBySectors();
    }

    double measure (double beta, const SMatrix& O)
    {
        return ensembleAverage(beta, O);
    }

    double measurePolarization2 (double beta, size_t p)
    {
        /*
         *        (n_0+n_2)^2 - (n_0-n_2)^2
         * d_02 = -------------------------
         *               (n_0+n_2)^2
         *
         * d_02 measures how evenly distributed the charges are along the
         * diagonal 02. d_02 = 1 => evenly distributed. d_02 = 0 => unevenly
         * distributed.
         *
         * P = d_02 * d_13 * 1/2 * (n_1 + n_3 - n_0 - n_2)
         */
        const size_t o = 4*p;
        const double n0 = ensembleAverage(beta, s.n(o+0));
        const double n1 = ensembleAverage(beta, s.n(o+1));
        const double n2 = ensembleAverage(beta, s.n(o+2));
        const double n3 = ensembleAverage(beta, s.n(o+3));
        return 8 * n0*n1*n2*n3 * (n1+n3-n0-n2) / 
               ( (n0+n2)*(n0+n2) * (n1+n3)*(n1+n3) );
    }

    std::vector<double> measureParticleNumber (double beta, size_t p)
    {
        const size_t o = 4*p;
        const double n0 = ensembleAverage(beta, s.n(o+0));
        const double n1 = ensembleAverage(beta, s.n(o+1));
        const double n2 = ensembleAverage(beta, s.n(o+2));
        const double n3 = ensembleAverage(beta, s.n(o+3));
        std::vector<double> ns;
        ns.push_back(n0);
        ns.push_back(n1);
        ns.push_back(n2);
        ns.push_back(n3);
        ns.push_back(n0+n1+n2+n3);
        return ns;
    }

    const DVector& energies ()
    {
        return H.eigenvalues();
    }

    double Emin () const
    {
        return H.Emin();
    }

    size_t N_p, N_sites, electronsPerPlaquet;
    QcaHamiltonian<S> H;
    EnsembleAverageBySectors<S> ensembleAverage;
    Polarization<S> P;
    ParticleNumber<S> N;
    double t, td, ti, V0, a, b, Vext, Pext, mu, epsilonr, lambdaD, epsilon0, q;
};

template<template <typename> class ExternalTC>
class QcaBond : public QcaCommon< QcaBond<ExternalTC> >
{
public:
    typedef QcaBond<ExternalTC> Self;
    typedef QcaCommon<Self> Base;
    typedef ExternalTC<Self> External;

    QcaBond (size_t N_p_)
        : Base(*this, N_p_), basis(plaquetSize*N_p_), ca(*this, plaquetSize), 
          PPSO(plaquetSize)
    {
        basis.addSymmetryOperator(&PPSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(2,Base::N_p);
        basis.setFilter(constructSector(filterValue));
        basis.construct();
        ca.construct();
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
    CreatorAnnihilator<QcaBond> ca;

private:
    ParticleNumberPerPlaquetSymmetryOperator PPSO;
};

template<template <typename> class ExternalTC, size_t numberOfElectronsPerPlaquet = 2>
class QcaFixedCharge : public QcaCommon< QcaFixedCharge<ExternalTC, numberOfElectronsPerPlaquet> >
{
public:
    typedef QcaFixedCharge<ExternalTC, numberOfElectronsPerPlaquet> Self;
    typedef QcaCommon<Self> Base;
    typedef ExternalTC<Self> External;

    QcaFixedCharge (size_t N_p_)
        : Base(*this, N_p_, numberOfElectronsPerPlaquet), basis(plaquetSize*N_p_), 
          creatorAnnihilator(*this, plaquetSize), PPSO(plaquetSize)
    {
        basis.addSymmetryOperator(&PPSO);
        basis.addSymmetryOperator(&SSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(numberOfElectronsPerPlaquet, Base::N_p);
        basis.setFilter(constructSector(filterValue));
        basis.construct();
        creatorAnnihilator.construct();
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
    CreatorAnnihilator<QcaFixedCharge> creatorAnnihilator;

private:
    ParticleNumberPerPlaquetSymmetryOperator PPSO;
    SpinSymmetryOperator SSO;
};

template<template <typename> class ExternalTC, size_t numberOfElectronsPerPlaquet = 2>
class QcaGrandCanonical : public QcaCommon< QcaGrandCanonical<ExternalTC, numberOfElectronsPerPlaquet> >
{
public:
    typedef QcaGrandCanonical<ExternalTC, numberOfElectronsPerPlaquet> Self;
    typedef QcaCommon<Self> Base;
    typedef ExternalTC<Self> External;

    QcaGrandCanonical (size_t N_p_)
        : Base(*this, N_p_, numberOfElectronsPerPlaquet), basis(plaquetSize*N_p_), creator(*this), 
          annihilator(*this) 
    {
        basis.addSymmetryOperator(&PSO);
        basis.addSymmetryOperator(&SSO);
        basis.construct();
        creator.construct();
        annihilator.construct();
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
    Creator<QcaGrandCanonical> creator;
    Annihilator<QcaGrandCanonical> annihilator;

private:
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
        QcaSystem::t = os["t"].get<double>(1.0);
        QcaSystem::td = os["td"].get<double>(0); 
        QcaSystem::ti = os["ti"].get<double>(0); 
        QcaSystem::a = os["a"].get<double>(1.0); 
        QcaSystem::b = os["b"].get<double>(3);
        QcaSystem::Vext = os["Vext"].get<double>(0);
        QcaSystem::Pext = os["Pext"].get<double>(0);
        QcaSystem::V0 = os["V0"].get<double>(1000); 
        QcaSystem::mu = os["mu"].get<double>(0);
        QcaSystem::epsilonr = os["epsilonr"].get<double>(1);
        QcaSystem::lambdaD = os["lambdaD"].get<double>(0);
        QcaSystem::q = os["q"].get<double>(0);
    }

    OptionSection getParameters ()
    {
        OptionSection os;
        os["p"] = QcaSystem::N_p;
        os["t"] = QcaSystem::t;
        os["td"] = QcaSystem::td;
        os["ti"] = QcaSystem::ti;
        os["V0"] = QcaSystem::V0;
        os["mu"] = QcaSystem::mu;
        os["Vext"] = QcaSystem::Vext;
        os["Pext"] = QcaSystem::Pext;
        os["a"] = QcaSystem::a;
        os["b"] = QcaSystem::b;
        os["epsilonr"] = QcaSystem::epsilonr;
        os["lambdaD"] = QcaSystem::lambdaD;
        os["q"] = QcaSystem::q;

        return os;
    }
};

/*
 * Useful typedefs
 */
typedef QcaBond<ExternalPlain> QcaBondPlain;
typedef QcaBond<ExternalDeadPlaquet> QcaBondDeadPlaquet;
typedef QcaFixedCharge<ExternalPlain, 2> QcaFixedCharge2Plain;
typedef QcaFixedCharge<ExternalPlain, 6> QcaFixedCharge6Plain;
typedef QcaFixedCharge<ExternalDeadPlaquet, 2> QcaFixedCharge2DeadPlaquet;
typedef QcaFixedCharge<ExternalDeadPlaquet, 6> QcaFixedCharge6DeadPlaquet;
typedef QcaGrandCanonical<ExternalPlain> QcaGrandCanonicalPlain;
typedef QcaGrandCanonical<ExternalDeadPlaquet, 2> QcaGrandCanonical2DeadPlaquet;
typedef QcaGrandCanonical<ExternalDeadPlaquet, 6> QcaGrandCanonical6DeadPlaquet;

typedef DQcaGeneric<QcaBond<ExternalPlain> > DQcaBondPlain;
typedef DQcaGeneric<QcaBond<ExternalDeadPlaquet> > DQcaBondDeadPlaquet;
typedef DQcaGeneric<QcaFixedCharge<ExternalPlain, 2> > DQcaFixedCharge2Plain;
typedef DQcaGeneric<QcaFixedCharge<ExternalPlain, 6> > DQcaFixedCharge6Plain;
typedef DQcaGeneric<QcaFixedCharge<ExternalDeadPlaquet, 2> > DQcaFixedCharge2DeadPlaquet;
typedef DQcaGeneric<QcaFixedCharge<ExternalDeadPlaquet, 6> > DQcaFixedCharge6DeadPlaquet;
typedef DQcaGeneric<QcaGrandCanonical<ExternalPlain> > DQcaGrandCanonicalPlain;
typedef DQcaGeneric<QcaGrandCanonical<ExternalDeadPlaquet, 2> > DQcaGrandCanonical2DeadPlaquet;
typedef DQcaGeneric<QcaGrandCanonical<ExternalDeadPlaquet, 6> > DQcaGrandCanonical6DeadPlaquet;

#endif // __QCA_HPP__
