#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"
#include "utilities.hpp"
#include <limits>

const double QCA_ELEMENTARY_CHARGE = 1.602176565E-19;
const double QCA_EPSILON_0 = 8.8541878176E-12;
const double QCA_NATURAL_EPSILON_R = QCA_ELEMENTARY_CHARGE / (4*M_PI*QCA_EPSILON_0*1e-9);
enum ElectronsPerCell {epc2 = 2, epc6 = 6};

// TODO: for the header-part of the printout -- how to best report the layout
// used; this includes: a, b, Pext
class Layout
{
private:
    std::vector<Vector2d> r_sites;
    std::vector<Vector2d> r_charges;
    std::vector<double> charges;
    ElectronsPerCell epc;
public:
    Layout (ElectronsPerCell epc_ = epc2)
    : epc(epc_)
    {}

    Layout& addSite (double r_x, double r_y)
    {
        r_sites.push_back(Vector2d(r_x, r_y));
        return *this;
    }

    Layout& addCharge (double r_x, double r_y, double c)
    {
        r_charges.push_back(Vector2d(r_x, r_y));
        charges.push_back(c);
        return *this;
    }

    Layout& addCell (double r_x, double r_y, double a)
    {
        addSite(r_x, r_y);
        addSite(r_x, r_y+a);
        addSite(r_x+a, r_y+a);
        addSite(r_x+a, r_y);
        return *this;
    }

    Layout& addDriverCell (double r_x, double r_y, double a, double P)
    {
        // compensation charge. 
        // q=0 for 2 electrons per cell, q=1 for 6 electrons per cell
        assert(epc==2 || epc==6);
        double q = 0;
        if (epc==2) q=0;
        if (epc==6) q=1;
        addCharge(r_x, r_y, q + (P+1)/2);
        addCharge(r_x, r_y+a, q + (1-P)/2);
        addCharge(r_x+a, r_y+a, q + (P+1)/2);
        addCharge(r_x+a, r_y, q + (1-P)/2);
        return *this;
    }

    Layout& addWire (double r_x, double r_y, int N_p, double a, double b, double P)
    {
        addDriverCell (r_x-b-a, r_y, a, P);
        for (int i=0; i<N_p; i++)
            addCell(r_x+i*(a+b), r_y, a);
        return *this;
    }

    Layout& addWireNoDriver (double r_x, double r_y, int N_p, double a, double b)
    {
        for (int i=0; i<N_p; i++)
            addCell(r_x+i*(a+b), r_y, a);
        return *this;
    }

    Layout& addNonuniformWire (double r_x, double r_y, int N_p, double a, 
                               std::vector<double> bs, double P)
    {
        assert (N_p == static_cast<int>(bs.size()));
        addDriverCell (r_x-bs[0]-a, r_y, a, P);
        double x_off=0;
        for (int i=0; i<static_cast<int>(bs.size()); i++)
        {
            if (i!=0) x_off += a+bs[i];
            addCell(x_off+r_x, r_y, a);
        }
        return *this;
    }

    void wire (int N_p, double a, double b, double P, ElectronsPerCell epc_)
    {
        clear();
        setElectronsPerCell(epc_);
        addWire(0,0, N_p, a, b, P);
    }

    void wire2e (int N_p, double a=1, double b=3, double P=0)
    {
        wire(N_p, a, b, P, epc2);
    }

    void wire6e (int N_p, double a=1, double b=3, double P=0)
    {
        wire(N_p, a, b, P, epc6);
    }

    void wireNoDriver (int N_p, double a, double b, ElectronsPerCell epc_)
    {
        clear();
        setElectronsPerCell(epc_);
        addWireNoDriver(0,0, N_p, a, b);
    }

    void wireNoDriver2e (int N_p, double a=1, double b=3)
    {
        wireNoDriver(N_p, a, b, epc2);
    }

    void wireNoDriver6e (int N_p, double a=1, double b=3)
    {
        wireNoDriver(N_p, a, b, epc6);
    }

    void nonuniformWire (int N_p, double a, std::vector<double> bs, double P, ElectronsPerCell epc_)
    {
        clear();
        setElectronsPerCell(epc_);
        addNonuniformWire(0,0, N_p, a, bs, P);
    }

    void nonuniformWire2e (int N_p, double a, std::vector<double> bs, double P)
    {
        nonuniformWire(N_p, a, bs, P, epc2);
    }

    void nonuniformWire6e (int N_p, double a, std::vector<double> bs, double P)
    {
        nonuniformWire(N_p, a, bs, P, epc6);
    }

    int N_sites () const
    {
        return static_cast<int>(r_sites.size());
    }

    int N_charges () const
    {
        assert (r_charges.size() == charges.size());
        return static_cast<int>(r_charges.size());
    }

    double r (int i, int j) const
    {
        const Vector2d d = r_sites[i] - r_sites[j];
        return d.norm();
    }

    double r_charge_dot (int i, int j) const
    {
        const Vector2d d = r_charges[i] - r_sites[j];
        return d.norm();
    }

    double charge (int i) const
    {
        return charges[i];
    }

    int getElectronsPerCell () const
    {
        return epc;
    }

    void setElectronsPerCell (ElectronsPerCell epc_)
    {
        epc = epc_;
    }

    void clear ()
    {
        r_sites.clear();
        r_charges.clear();
        charges.clear();
    }
};

class Wire2e : public Layout
{
public:
    Wire2e (int N_p, double a = 1, double b = 3, double P = 0)
    {
        wire2e(N_p, a, b, P);
    }
};

class Wire6e : public Layout
{
public:
    Wire6e (int N_p, double a = 1, double b = 3, double P = 0)
    {
        wire6e(N_p, a, b, P);
    }
};

class WireNoDriver2e : public Layout
{
public:
    WireNoDriver2e (int N_p, double a = 1, double b = 3)
    {
        wireNoDriver2e(N_p, a, b);
    }
};

class WireNoDriver6e : public Layout
{
public:
    WireNoDriver6e (int N_p, double a = 1, double b = 3)
    {
        wireNoDriver6e(N_p, a, b);
    }
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
        constructIdentityMatrix();
        H = SMatrix(s.basis.size(), s.basis.size());
        H.setZero();
        for (int i=0; i<s.N_sites(); i++)
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
            //TODO: What about the compensation charge? Should it not be n-q for
            //the external term -- this might be a serious bug
            H += (external(i) - s.mu) * s.n(i);
            for (int j=i+1; j<s.N_sites(); j++)
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
    
    //TODO: rewrite this so that hopping depends on the distance, i.e. t_ij =
    //t_ij(r_ij)
    double hopping (int i, int j) const
    {
        /*
         * Same cell
         */
        if (i/4 == j/4)
        {
            if (std::abs(i-j) == 2)
                return s.td;
            else if (i != j)
                return s.t;
        }
        /*
         * No inter-cell hopping for now
         */
        return 0;
    }
    
    double coulomb (int i, int j) const
    {
        if (i == j)
            return s.V0;
        const double r = s.l.r(i,j);
        if (s.lambdaD == 0)
            return QCA_ELEMENTARY_CHARGE / 
                   (4*M_PI * s.epsilon0 * s.epsilonr * r * 1e-9);
        else
            return QCA_ELEMENTARY_CHARGE * exp(- r / s.lambdaD) / 
                   (4*M_PI * s.epsilon0 * s.epsilonr * r * 1e-9);
    }

    double external (int i) const
    {
        double V=0;
        
        // simple external potential
        if (i==1)
            V += s.Vext;
        if (i==0)
            V += -s.Vext;

        // external potential due to static charges, e.g. a driver cell
        // formerly I called this dead plaquet
        for (int j=0; j<s.l.N_charges(); j++)
        {
            const double r = s.l.r_charge_dot(j,i);
            if (s.lambdaD == 0)
                V += (s.l.charge(j) - s.q) * QCA_ELEMENTARY_CHARGE / 
                       (4*M_PI * s.epsilon0 * s.epsilonr * r * 1e-9);
            else
                V += (s.l.charge(j) - s.q) * QCA_ELEMENTARY_CHARGE * exp(- r / s.lambdaD) / 
                       (4*M_PI * s.epsilon0 * s.epsilonr * r * 1e-9);
        }
        
        return V;
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
        return 1.0/2.0 * ( s.n(o+0)+s.n(o+2) - s.n(o+1)-s.n(o+3) );
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
        for (int i=0; i<s.N_sites(); i++)
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
    CreatorAnnihilator (const System& s_)
    : s(s_)
    {}

    void construct ()
    {
        cas = std::vector<SMatrix>(s.N_p() * s.plaquetSize * s.plaquetSize);
        zeroMatrix = SMatrix(s.basis.size(), s.basis.size());
        //TODO: optimise - c_i a_j = (c_j a_i)^{\dag}
        for (int p=0; p<s.N_p(); p++)
            for (int i=0; i<s.plaquetSize; i++)
                for (int j=0; j<s.plaquetSize; j++)
                    constructMatrix(s.plaquetSize*p + i, s.plaquetSize*p + j);
    }

    const SMatrix& operator() (size_t i, size_t j) const
    {
        // no inter-cell hopping => return a 0-matrix
        if (i/s.plaquetSize != j/s.plaquetSize)
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
        const size_t p = i/s.plaquetSize;
        assert(p == j/s.plaquetSize);
        const size_t ii = i%s.plaquetSize;
        const size_t jj = j%s.plaquetSize;
        return s.plaquetSize * s.plaquetSize * p + s.plaquetSize * ii + jj;
    }

    const System& s;
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

    int valueForNElectronsPerPlaquet (int N, int N_p) const
    {
        assert(N_p <= std::numeric_limits<int>::digits10);
        int value = 0;
        int multiplier = 1;
        for (int i=0; i<N_p; i++)
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
    int N_p_, N_sites_;

public:
    QcaHamiltonian<S> H;
    EnsembleAverageBySectors<S> ensembleAverage;
    Polarization<S> P;
    ParticleNumber<S> N;
    Layout l;
    double t, td, V0, Vext, mu, epsilonr, lambdaD, epsilon0, q, beta;

public:
    QcaCommon (QcaSystem& s_)
        : s(s_), N_p_(0), N_sites_(0), 
          H(s), ensembleAverage(s), P(s), N(s), 
          t(1), td(0), V0(1000), Vext(0), mu(0),
          epsilonr(QCA_NATURAL_EPSILON_R), lambdaD(0), 
          epsilon0(QCA_EPSILON_0), q(0), beta(1)
    {}

    void update ()
    {
        if (l.N_sites() != N_sites_)
            s.constructBasis();
        H.construct();
        H.diagonalizeUsingSymmetriesBySectors();
    }

    double measure (double beta_, const SMatrix& O) const
    {
        return ensembleAverage(beta_, O);
    }

    double measure (const SMatrix& O) const
    {
        return measure(beta, O);
    }

    double measurePolarization (double beta_, int p) const
    {
        return ensembleAverage(beta_, P(p));
    }

    double measurePolarization (int p) const
    {
        return measurePolarization(beta, p);
    }

    double measurePolarization2 (double beta_, int p) const
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
        const int o = 4*p;
        const double n0 = ensembleAverage(beta_, s.n(o+0));
        const double n1 = ensembleAverage(beta_, s.n(o+1));
        const double n2 = ensembleAverage(beta_, s.n(o+2));
        const double n3 = ensembleAverage(beta_, s.n(o+3));
        return 8 * n0*n1*n2*n3 * (n0+n2-n1-n3) / 
               ( (n0+n2)*(n0+n2) * (n1+n3)*(n1+n3) );
    }
    
    double measurePolarization2 (int p) const
    {
        return measurePolarization2(beta, p);
    }

    std::vector<double> measureParticleNumber (double beta_, int p) const
    {
        const int o = 4*p;
        const double n0 = ensembleAverage(beta_, s.n(o+0));
        const double n1 = ensembleAverage(beta_, s.n(o+1));
        const double n2 = ensembleAverage(beta_, s.n(o+2));
        const double n3 = ensembleAverage(beta_, s.n(o+3));
        std::vector<double> ns;
        ns.push_back(n0);
        ns.push_back(n1);
        ns.push_back(n2);
        ns.push_back(n3);
        ns.push_back(n0+n1+n2+n3);
        return ns;
    }
    
    std::vector<double> measureParticleNumber (int p) const
    {
        return measureParticleNumber(beta, p);
    }

    //TODO: better return a std::vector -- only use Eigen internally (inside Qca
    //classes)
    const DVector& energies ()
    {
        return H.eigenvalues();
    }

    double Emin () const
    {
        return H.Emin();
    }

    int N_p () const
    {
        return N_p_;
    }

    int N_sites () const
    {
        return N_sites_;
    }

protected:
    void updateParametersFromLayout ()
    {
            N_sites_ = l.N_sites();
            N_p_ = l.N_sites()/4;
            assert(N_sites_ = N_p_ * 4);
            assert(l.N_charges() == 4 || l.N_charges() == 0);
    }
};

class QcaBond : public QcaCommon<QcaBond>
{
public:
    typedef QcaBond Self;
    typedef QcaCommon<Self> Base;

    enum {plaquetSize=4};
    Basis basis;
    CreatorAnnihilator<Self> ca;

private:
    ParticleNumberPerPlaquetSymmetryOperator PPSO;

public:
    QcaBond (Layout l_ = Layout())
        : Base(*this), ca(*this), PPSO(plaquetSize)
    {
        Base::l = l_;
    }

    void constructBasis ()
    {
        Base::updateParametersFromLayout();
        basis = Basis();
        basis.addSymmetryOperator(&PPSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(2,Base::N_p());
        basis.setFilter(constructSector(filterValue));
        basis.construct(plaquetSize*Base::N_p());
        ca.construct();
    }

    SMatrix n (int i) const
    {
        return ca(i,i);
    }

    SMatrix n_updown (int i) const
    {
        // return 0
        return SMatrix(basis.size(), basis.size());
    }
};

class QcaFixedCharge : public QcaCommon<QcaFixedCharge>
{
public:
    typedef QcaFixedCharge Self;
    typedef QcaCommon<Self> Base;

    enum {plaquetSize=8};
    Basis basis;
    CreatorAnnihilator<Self> creatorAnnihilator;

private:
    ParticleNumberPerPlaquetSymmetryOperator PPSO;
    SpinSymmetryOperator SSO;

public:
    QcaFixedCharge (Layout l_ = Layout())
        : Base(*this), creatorAnnihilator(*this), PPSO(plaquetSize)
    {
        Base::l = l_;
    }

    void constructBasis ()
    {
        Base::updateParametersFromLayout();
        basis = Basis();
        basis.addSymmetryOperator(&PPSO);
        basis.addSymmetryOperator(&SSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(l.getElectronsPerCell(), Base::N_p());
        basis.setFilter(constructSector(filterValue));
        basis.construct(plaquetSize*Base::N_p());
        creatorAnnihilator.construct();
    }

    size_t I (int i, Spin s) const
    {
        return 2*i + s;
    }

    SMatrix ca (int i, Spin s_i, int j, Spin s_j) const
    {
        return creatorAnnihilator(I(i, s_i), I(j, s_j));
    }

    SMatrix ca (int i, int j) const
    {
        return ca(i,UP,j,UP) + ca(i,DOWN,j,DOWN);
    }

    SMatrix n (int i, Spin s) const
    {
        return ca(i,s,i,s);
    }

    SMatrix n (int i) const
    {
        return n(i,UP) + n(i,DOWN);
    }

    SMatrix n_updown (int i) const
    {
        return n(i,UP) * n(i,DOWN);
    }
};

class QcaGrandCanonical : public QcaCommon<QcaGrandCanonical>
{
public:
    typedef QcaGrandCanonical Self;
    typedef QcaCommon<Self> Base;

    enum {plaquetSize=8};
    Basis basis;
    Creator<Self> creator;
    Annihilator<Self> annihilator;

private:
    ParticleNumberSymmetryOperator PSO;
    SpinSymmetryOperator SSO;

public:
    QcaGrandCanonical (Layout l_ = Layout())
        : Base(*this), creator(*this), annihilator(*this) 
    {
        Base::l = l_;
    }

    void constructBasis ()
    {
        Base::updateParametersFromLayout();
        basis = Basis();
        basis.addSymmetryOperator(&PSO);
        basis.addSymmetryOperator(&SSO);
        basis.construct(plaquetSize*Base::N_p());
        creator.construct();
        annihilator.construct();
    }

    size_t I (int i, Spin s) const
    {
        return 2*i + s;
    }

    SMatrix ca (int i, Spin s_i, int j, Spin s_j) const
    {
        return creator(I(i, s_i))*annihilator(I(j, s_j));
    }

    SMatrix ca (int i, int j) const
    {
        return ca(i,UP,j,UP) + ca(i,DOWN,j,DOWN);
    }

    SMatrix n (int i, Spin s) const
    {
        return ca(i,s,i,s);
    }

    SMatrix n (int i) const
    {
        return n(i,UP) + n(i,DOWN);
    }

    SMatrix n_updown (int i) const
    {
        return n(i,UP) * n(i,DOWN);
    }
};

template<class QcaSystem>
class DQcaGeneric : public QcaSystem
{
private:
    typedef QcaSystem Base;

public:
    DQcaGeneric (OptionSection os)
    : QcaSystem()
    {
        setParameters(os);
    }

    void setParameters (OptionSection os)
    {
        Base::t = os["t"].get<double>(1.0);
        Base::td = os["td"].get<double>(0); 
        Base::Vext = os["Vext"].get<double>(0);
        Base::V0 = os["V0"].get<double>(1000); 
        Base::mu = os["mu"].get<double>(0);
        Base::epsilonr = os["epsilonr"].get<double>(1);
        Base::lambdaD = os["lambdaD"].get<double>(0);
        Base::q = os["q"].get<double>(0);
        
        ElectronsPerCell epc = epc2;
        if (os["epc"].get<int>() == 2)
            epc = epc2;
        else if (os["epc"].get<int>() == 6)
            epc = epc6;
        
        if (os["layout"] == "wire")
            Base::l.wire(os["p"].get<int>(1), 
                         os["a"].get<double>(1.0), 
                         os["b"].get<double>(3.0), 
                         os["Pext"].get<double>(0), 
                         epc);
        else if (os["layout"] == "nonuniformwire")
        {
            std::vector<double> bs = os["bs"];
            Base::l.nonuniformWire(os["p"].get<int>(1), 
                                   os["a"].get<double>(1.0), 
                                   bs,
                                   os["Pext"].get<double>(0), 
                                   epc);
        }
        // TODO: better error handling
    }

    OptionSection getParameters ()
    {
        OptionSection os;
        os["p"] = Base::N_p();
        os["t"] = Base::t;
        os["td"] = Base::td;
        os["V0"] = Base::V0;
        os["mu"] = Base::mu;
        os["Vext"] = Base::Vext;
        os["epsilonr"] = Base::epsilonr;
        os["lambdaD"] = Base::lambdaD;
        os["q"] = Base::q;
        //TODO: read back layout parameters or description
        //os["a"] = Base::a;
        //os["b"] = Base::b;
        //os["Pext"] = Base::Pext;

        return os;
    }
};

#endif // __QCA_HPP__
