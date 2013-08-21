#ifndef __QCA_HPP__
#define __QCA_HPP__

#include "system.hpp"
#include "utilities.hpp"
#include <limits>

const double QCA_ELEMENTARY_CHARGE = 1.602176565E-19;
const double QCA_EPSILON_0 = 8.8541878176E-12;
const double QCA_NATURAL_EPSILON_R = QCA_ELEMENTARY_CHARGE / (4*M_PI*QCA_EPSILON_0*1e-9);
enum ElectronsPerCell {epc2 = 2, epc6 = 6};

class Layout
{
public:
    std::vector<Vector2d> r_sites;
    std::vector<Vector2d> r_charges;
    std::vector<double> charges;
    ElectronsPerCell epc;
public:
    Layout ()
    : epc(epc2)
    {}

    void addSite (double r_x, double r_y)
    {
        r_sites.push_back(Vector2d(r_x, r_y));
    }

    void addCharge (double r_x, double r_y, double c)
    {
        r_charges.push_back(Vector2d(r_x, r_y));
        charges.push_back(c);
    }

    void addCell (double r_x, double r_y, double a)
    {
        addSite(r_x, r_y);
        addSite(r_x, r_y+a);
        addSite(r_x+a, r_y+a);
        addSite(r_x+a, r_y);
    }

    void addDriverCell (double r_x, double r_y, double a, double P, ElectronsPerCell epc_)
    {
        // compensation charge. 
        // q=0 for 2 electrons per cell, q=1 for 6 electrons per cell
        assert(epc_==2 || epc_==6);
        double q = 0;
        if (epc_==2) q=0;
        if (epc_==6) q=1;
        addCharge(r_x, r_y, q + (P+1)/2);
        addCharge(r_x, r_y+a, q + (1-P)/2);
        addCharge(r_x+a, r_y+a, q + (P+1)/2);
        addCharge(r_x+a, r_y, q + (1-P)/2);
    }

    void addDriverCell (double r_x, double r_y, double a, double P)
    {
        return addDriverCell(r_x, r_y, a, P, epc);
    }

    void wire (int N_p, double a, double b, double P, ElectronsPerCell epc_)
    {
        clear();

        double r_x = 0;
        double r_y = 0;

        addDriverCell(r_x-b-a, r_y, a, P, epc_);
        for (int i=0; i<N_p; i++)
            addCell(r_x+i*(a+b), r_y, a);
    }

    void wire (int N_p, double a, double b, double P)
    {
        wire(N_p, a, b, P, epc);
    }

    void nonuniformWire (int N_p, double a, std::vector<double> bs, double P, ElectronsPerCell epc_)
    {
        clear();

        double r_x = 0;
        double r_y = 0;

        assert (N_p == static_cast<int>(bs.size()));
        addDriverCell(r_x-bs[0]-a, r_y, a, P, epc_);
        double x_off=0;
        for (int i=0; i<static_cast<int>(bs.size()); i++)
        {
            if (i!=0) x_off += a+bs[i];
            addCell(x_off+r_x, r_y, a);
        }
    }

    void nonuniformWire (int N_p, double a, std::vector<double> bs, double P)
    {
        nonuniformWire(N_p, a, bs, P, epc);
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

    void clear ()
    {
        r_sites.clear();
        r_charges.clear();
        charges.clear();
    }

    bool operator== (const Layout& l) const
    {
        return 
            r_sites == l.r_sites &&
            r_charges == l.r_charges &&
            charges == l.charges;
    }
};

template<class System>
class QcaHamiltonian : public Hamiltonian<System>
{
public:
    SMatrix& H;
    const System& s;

    QcaHamiltonian (const System& s_)
    : Hamiltonian<System>(s_), H(Hamiltonian<System>::H), 
      s(Hamiltonian<System>::s)
    {}

    void construct() 
    {
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
                H += coulomb(i,j) * ( s.n(i) * s.n(j) - s.q * ( s.n(i) + s.n(j) ) );
            }
        }
    }

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
        // external potential due to static charges, e.g. a driver cell
        // formerly I called this dead plaquet
        double V=0;
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
        m.reserve(VectorXi::Constant(s.basis.size(), 1));
        for (size_t col=0; col<s.basis.size(); col++)
        {
            if (i==j && s.basis(col)[i] == 1)
            {
                m.insert(col, col) = 1;
                continue;
            }
            if (s.basis(col)[i] == 1 || s.basis(col)[j] == 0)
                continue;
            State state(s.basis(col));
            state[i] = 1;
            state[j] = 0;
            const size_t row = s.basis(state);
            size_t sum = state.count(i,j);
            if (i<j) sum -= 1; //works, because for i<j we always have sum>=1
            const double sign = (sum%2==0)?1:-1; // probably faster than using (-1)^sum
            m.insert(row, col) = sign;
        }
        m.makeCompressed();
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

template<class QcaSystem, template <typename System> class Hamiltonian_=QcaHamiltonian>
class QcaCommon
{
private:
    typedef QcaSystem S;
    S& s;
    int N_p_, N_sites_;

public:
    Hamiltonian_<S> H;
    EnsembleAverageBySectors<S> ensembleAverage;
    Polarization<S> P;
    ParticleNumber<S> N;
    Layout l;
    double t, td, V0, mu, epsilonr, lambdaD, epsilon0, q, beta;

public:
    QcaCommon (QcaSystem& s_)
        : s(s_), N_p_(0), N_sites_(0), 
          H(s), ensembleAverage(s), P(s), N(s), l(Layout()),
          t(1), td(0), V0(1000), mu(0),
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

    std::vector<std::vector<double>> measureParticleNumberOverEnergy ()
    {
        if (s.N_sites()==0)
            return std::vector<std::vector<double>>();

        // Construct the operator O which measures the overall particle number
        SMatrix O = s.n(0);
        for (int i=1; i<s.N_sites(); i++)
            O += s.n(i);

        // Partition function
        double Z = 0;
        const std::vector<DVector>& eigenvalues = H.eigenvaluesBySectors();
        const std::vector<DMatrix>& eigenvectors = H.eigenvectorsBySectors();
        for (size_t i=0; i<eigenvalues.size(); i++)
            for (int j=0; j<eigenvalues[i].size(); j++)
                Z += std::exp(-beta * (eigenvalues[i](j) - H.Emin()));

        // Calculate the particle number / occupancy of each energy level
        std::vector<std::vector<double>> Ns(H.eigenvalues().size(), std::vector<double>(2));
        size_t index = 0;
        for (size_t i=0; i<eigenvalues.size(); i++)
        {
            const int size = eigenvalues[i].size();
            SMatrix O_block = EigenHelpers::sparseToSparseBlock(O, index, index, size, size);
            for (int j=0; j<size; j++)
            {
                Ns[index+j][0] = eigenvalues[i](j);
                Ns[index+j][1] = 
                    std::exp(-beta * (eigenvalues[i](j) - H.Emin())) / Z * 
                    eigenvectors[i].col(j).adjoint() * O_block * eigenvectors[i].col(j);
            }
            index += size;
        }
        return Ns;
    }

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
    QcaBond ()
        : Base(*this), ca(*this), PPSO(plaquetSize)
    {}

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
    QcaFixedCharge ()
        : Base(*this), creatorAnnihilator(*this), PPSO(plaquetSize)
    {}

    void constructBasis ()
    {
        Base::updateParametersFromLayout();
        basis = Basis();
        basis.addSymmetryOperator(&PPSO);
        basis.addSymmetryOperator(&SSO);
        int filterValue = PPSO.valueForNElectronsPerPlaquet(Base::l.epc, Base::N_p());
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
    QcaGrandCanonical ()
        : Base(*this), creator(*this), annihilator(*this) 
    {}

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

template<class System>
class QcaIsingHamiltonian : public QcaHamiltonian<System>
{
public:
    typedef QcaHamiltonian<System> Base;
    SMatrix& H;
    const System& s;

    QcaIsingHamiltonian (const System& s_)
    : Base(s_), H(Base::H), s(Base::s)
    {}

    void construct() 
    {
        /*
         * Construct the hopping part first.
         */
        H = SMatrix(s.basis.size(), s.basis.size());
        // We expect one entry per column
        H.reserve(VectorXi::Constant(s.basis.size(), 1));
        for (size_t col=0; col<s.basis.size(); col++)
        {
            State state(s.basis(col));
            for (size_t i=0; i<s.N_p(); i++)
            {
                // flip spin: 1 -> 0, 0 -> 1
                state[i] = 1 - state[i];
            }
            const size_t row = s.basis(state);
            // Note: this t is different from the ts in the other QCA
            // Hamiltonians!
            H.insert(row,col) = - s.t * s.N_p();
        }
        H.makeCompressed();

        /*
         * Add external potential and Coulomb interaction.
         */
        for (int i=0; i<s.N_sites(); i++)
        {
            /*
             * Eigen seems to truncate very small values in Sparse matrices. The
             * epsilon value used for the truncation is defined in Eigen's
             * NumTraits. To be safe we check that our values are bigger than
             * this threshold. 
             */
            assert(Base::coulomb(i,i)==0 || 
                   std::fabs(Base::coulomb(i,i))>NumTraits<double>::dummy_precision());
            assert(Base::external(i)==0 || 
                   std::fabs((Base::external(i)+s.mu))> NumTraits<double>::dummy_precision());

            H += Base::external(i) * s.n(i);
            for (int j=i+1; j<s.N_sites(); j++)
                H += Base::coulomb(i,j) * ( s.n(i) * s.n(j) - s.q * ( s.n(i) + s.n(j) ) );
        }
    }
};

template <class System>
class Sigma
{
private:
    const System& s;
    std::vector<SMatrix> ss;

public:
    Sigma (const System& s_)
    : s(s_)
    {}

    void construct()
    {
        ss.resize(s.N_p());
        for (size_t i=0; i<s.N_p(); i++)
            constructMatrix(i);
    }

    const SMatrix& operator() (size_t i) const
    {
        return ss[i];
    }

    /** Each matrix measures the "spin" (+1/-1) on cell/plaquet i. */
    void constructMatrix(size_t i)
    {
        SMatrix& m = ss[i];
        m = SMatrix(s.basis.size(), s.basis.size());
        // This is diagonal matrix, so one entry per column
        m.reserve(VectorXi::Constant(s.basis.size(), 1));
        for (size_t j=0; j<s.basis.size(); j++)
        {
            const State& state = s.basis(j);
            if (state[i] == 1)
                m.insert(j,j) = +1; //"spin" up
            else
                m.insert(j,j) = -1; //"spin" down
        }
        m.makeCompressed();
    }
};

class QcaIsing : public QcaCommon<QcaIsing, QcaIsingHamiltonian>
{
public:
    typedef QcaIsing Self;
    typedef QcaCommon<Self, QcaIsingHamiltonian> Base;
    Basis basis;
    Sigma<Self> sigma;

    QcaIsing()
    : Base(*this), sigma(*this)
    {}

    void constructBasis ()
    {
        Base::updateParametersFromLayout();
        basis = Basis();
        basis.construct(Base::N_p());
        sigma.construct();
    }

    SMatrix n (int i) const
    {
        // construct an identity matrix
        SMatrix I(basis.size(), basis.size());
        I.reserve(VectorXi::Constant(basis.size(),1));
        for (size_t k=0; k<basis.size(); k++)
            I.insert(k,k) = 1;
        I.makeCompressed();

        int c = i/4; // which cell
        int j = i%4; // which site on the cell
        if (j==0 || j==2)
            return 0.5 * (I + sigma(c));
        else // j==1 || j==3
            return 0.5 * (I - sigma(c));
    }

    double measureSpin (int i) const
    {
        return ensembleAverage(beta, sigma(i));
    }
};

#endif // __QCA_HPP__
