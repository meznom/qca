#ifndef __DIAGONALISER_HPP__
#define __DIAGONALISER_HPP__

#include <stdint.h>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <vector>
#include <cassert>
#include <cmath>
#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
using namespace Eigen;


//typedef std::vector<bool> state;
typedef std::bitset<64> state;

template<typename T=uint64_t>
class FermionicState
{
public:
    typedef T num_t; 
    typedef std::bitset<sizeof(T)*8> bitset;

    FermionicState (size_t N_) : N(N_)
    {
        assert(N<sizeof(T)*8);
    }

    void operator= (const T& v)
    {
        s = v;
    }

    bool operator== (const FermionicState<T>& fs) const
    {
        return s == fs.s && N == fs.N;
    }

    bool operator[] (size_t pos) const
    {
        return s[pos];
    }

    typename bitset::reference operator[] (size_t pos)
    {
        //return s[pos];
        return typename bitset::reference(s, pos);
    }

    size_t size() const
    {
        return N;
    }

    size_t count(size_t toPos) const
    {
        size_t sum = 0;
        for (size_t i=0; i<toPos; i++)
            sum += s[i];
        return sum;
    }

private:
    size_t N;
    bitset s;
};

template<typename T>
std::ostream& operator<< (std::ostream& o, const FermionicState<T>& fs)
{
    for (size_t i=0; i<fs.size(); i++)
        o << fs[i];
    return o;
}

class FilterSelectAll
{
public:
    bool operator() (const FermionicState<>&) {return true;}
};

class SorterDontSort
{
public:
    bool operator() (const FermionicState<>&, const FermionicState<>&) {return false;}
};

class System
{
public:
    typedef FermionicState<> State;
    typedef State::num_t num_t;
    typedef SparseMatrix<double> SMatrix;
    typedef MatrixXd DMatrix;
    typedef VectorXd DVector;
    typedef std::vector<State>::const_iterator CBIT;

    template<class Pred, class BinPred>
    System (size_t N_orbital_, Pred filter = FilterSelectAll(), BinPred sorter = SorterDontSort()) 
    : N_orbital(N_orbital_), cs(N_orbital_), as(N_orbital_), Emin(0)
    {
        //constructBasis(N_orbital, voidFilter, voidSorter);
        constructBasis(N_orbital, filter, sorter);

        H = SMatrix(basis.size(), basis.size());

        
        for (size_t i=0; i<N_orbital; i++)
        {
            cs[i] = creatorMatrix(i);
            as[i] = SMatrix(cs[i].transpose());
        }
    }

    template<class Pred, class BinPred>
    void constructBasis (size_t N_orbital, Pred filter, BinPred sorter)
    {
        N_basis = std::pow(2.0, static_cast<int>(N_orbital));
        basis = std::vector<State>(N_basis, State(N_orbital));
        State s(N_orbital);
        num_t i=0;
        for (num_t num=0; num<N_basis; num++)
        {
            s = num;
            if (filter(s))
            {
                basis[i] = num;
                i++;
            }
        }
        std::sort(basis.begin(), basis.end(), sorter);
    }

    State getStateFromIndex (num_t i) const
    {
        return basis[i];
    }

    num_t getIndexFromState (const State& s) const
    {
        for (CBIT i=basis.begin(); i!=basis.end(); i++)
        {
            if (*i == s)
                return i-basis.begin();
        }
        //TODO: throw exception if not found!
    }

    void dumpBasis ()
    {
        //for (CBIT i=basis.begin(); i!=basis.end(); i++)
        for (size_t i=0; i<basis.size(); i++)
        {
            std::cerr << basis[i] << std::endl;
        }
    }

    SMatrix creatorMatrix (size_t i)
    {
        SMatrix m(basis.size(), basis.size());
        // we expect one entry per column
        m.reserve(basis.size());
        for (size_t col=0; col<basis.size(); col++)
        {
            m.startVec(col);
            if (basis[col][i] == 1)
                continue;
            State s(basis[col]);
            s[i] = 1;
            const size_t row = getIndexFromState(s);
            const size_t sum = s.count(i);
            const double sign = (sum%2==0)?1:-1; // probably faster than using (-1)^sum
            m.insertBack(row, col) = sign;
        }
        m.finalize();

        return m;
    }

    size_t I (size_t i) const
    {
        return i;
    }

    const SMatrix& creator (size_t i) const
    {
        return cs[i];
    }

    const SMatrix& c (size_t i) const
    {
        return creator(I(i));
    }

    const SMatrix& annihilator (size_t i) const
    {
        return as[i];
    }

    const SMatrix& a (size_t i) const
    {
        return annihilator(I(i));
    }

    void constructH () {}

    void diagonalise ()
    {
        DMatrix m(H);
        SelfAdjointEigenSolver<DMatrix> es(m);
        eigenvalues = es.eigenvalues();
        eigenvectors = es.eigenvectors();
        Emin = eigenvalues.minCoeff();
    }

    double ensembleAverage (double beta, const SMatrix& O) const
    {
        double sum = 0;
        for (size_t i=0; i<eigenvalues.size(); i++)
            sum += 
                std::exp(-beta * (eigenvalues(i) - Emin)) * 
                eigenvectors.col(i).adjoint() * O * eigenvectors.col(i);
        return sum / partitionFunction(beta);
    }

    double partitionFunction (double beta) const
    {
        double Z = 0;
        for (size_t i=0; i<eigenvalues.size(); i++)
            Z += std::exp(-beta * (eigenvalues(i) - Emin));
        return Z;
    }

//private:
    std::vector<State> basis;
    num_t N_basis;
    size_t N_orbital;
    FilterSelectAll voidFilter;
    SorterDontSort voidSorter;
    std::vector<SMatrix> cs, as;
    enum Spin {UP=0, DOWN=1};
    
    SMatrix H;
    DMatrix eigenvectors;
    DVector eigenvalues;
    double Emin;
};

class AndersonModel : public System
{
public:
    AndersonModel (size_t N_bath_)
    : System(2*N_bath_+ 2, FilterSelectAll(), SorterDontSort()), N_bath(N_bath_), N_sites(N_bath_+1), t(1), mu(0), U(0), V(1), Ed(0)
    {}
    
    size_t I (size_t i, Spin s) const
    {
        return 2*i + s;
    }

    const SMatrix& c (size_t i, Spin s) const
    {
        return creator(I(i,s));
    }

    const SMatrix& a (size_t i, Spin s) const
    {
        return annihilator(I(i,s));
    }

    void constructH ()
    {
        H.setZero();
        for (size_t i = 1; i<N_sites; i++)
        {
            H += - mu * ( c(i,UP) * a(i,UP) + c(i,DOWN) * a(i,DOWN) );
            
            // this is physically correct, I think
            //size_t j = (i+1==N_sites)?1:i+1;
            //if (i==j) continue;
            //H += -t * ( 
            //    c(i,UP) * a(j,UP) + c(j,UP) * a(i,UP) +
            //    c(i,DOWN) * a(j,DOWN) + c(j,DOWN) * a(i,DOWN)
            //);

            // this is what my python diagonaliser does (and probably my ctqmc
            // code as well)
            for (size_t j=i-1; j<=i+1; j+=2)
            {
                size_t k = j;
                if (k==0) k = N_sites-1;
                else if (k==N_sites) k = 1;
                H += -t * ( 
                    c(i,UP) * a(k,UP) + c(i,DOWN) * a(k,DOWN)
                );
            }
        }
        H += Ed * (c(0,UP)*a(0,UP) + c(0,DOWN)*a(0,DOWN));
        H += V * (
            c(0,UP)*a(1,UP) + c(1,UP)*a(0,UP) + 
            c(0,DOWN)*a(1,DOWN) + c(1,DOWN)*a(0,DOWN)
        );
        //SMatrix I(DMatrix::Identity(basis.size(),basis.size()));
        //TODO
        SMatrix Id(basis.size(), basis.size());
        Id.reserve(basis.size());
        for (size_t i=0; i<basis.size(); i++)
        {
            Id.startVec(i);
            Id.insertBack(i,i) = 1;
        }
        Id.finalize();
        H += U * (c(0,UP)*a(0,UP) - 0.5*Id) * (c(0,DOWN)*a(0,DOWN) - 0.5*Id);
    }

    SMatrix doubleOccupancy (size_t i) const
    {
        return c(i,UP)*a(i,UP) * c(i,DOWN)*a(i,DOWN);
    }

    SMatrix n (size_t i, Spin s) const
    {
        return c(i,s)*a(i,s);
    }

    SMatrix n (size_t i) const
    {
        return n(i,UP) + n(i,DOWN);
    }

    SMatrix particleNumber () const
    {
        SMatrix pn(basis.size(), basis.size());
        for (size_t i=0; i<N_sites; i++)
        {
            pn += n(i);
        }
        return pn;
    }

    double t, mu, U, V, Ed;
    size_t N_bath, N_sites;
};

class Hopping
{
public:
    Hopping (double t_, double td_)
    : t(t_), td(td_)
    {}

    double operator() (size_t i, size_t j) const
    {
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
    : V0(0), a(a), b(b)
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
        assert( std::abs(deltaY==a) || deltaY==0 );

        /*
         * 0 1
         * 3 2
         */
        //TODO: check if this is correct, can we optimise it?
        const double deltaX = 
            (a+b) * ( (i - j) / 4 ) + 
            a * ( ( ((i%4)%3==0)?0:1 ) - ( ((j%4)%3==0)?0:1 ) );

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
    }
private:
    const double Vext;
};

class QCABond : public System
{
public:
    QCABond (size_t N_plaquet_)
    : System (4*N_plaquet_, FilterSelectAll(), SorterDontSort()), 
      N_plaquet(N_plaquet_), N_sites(4*N_plaquet_)
    {
        param.t = 1;
        param.td = 0;
        param.V0 = 1000;
        param.a = 1.0/200;
        param.b = 3*param.a;
        param.Vext = 0;
    }

    SMatrix ca (size_t i, size_t j) const
    {
        return c(i) * a(j);
    }

    SMatrix n (size_t i) const
    {
        return ca(i,i);
    }

    void constructH ()
    {
        Hopping hopping(param.t, param.td);
        Coulomb coulomb(param.V0, param.a, param.b);
        External external(param.Vext);

        H.setZero();
        for (size_t i=0; i<N_sites; i++)
        {
            H += external(i) * n(i);
            for (size_t j=i+1; j<N_sites; j++)
            {
                H += - hopping(i,j) * ca(i,j) - hopping(j,i) * ca(j,i);
                H += coulomb(i,j) * n(i) * n(j);
            }
        }
        
    }

    struct {
        double t, td, V0, a, b, Vext;
    } param;

private:
    size_t N_plaquet, N_sites;
};

template<class S>
class DoubleOccupancy
{
public:
    DoubleOccupancy (const S& s_) 
    : s(s_)
    {}

    typename S::SMatrix operator() (size_t i) const
    {
        return s.c(i,S::UP)*s.a(i,S::UP) * s.c(i,S::DOWN)*s.a(i,S::DOWN);
    }

private:
    const S& s;
};

template<class S>
class EnsembleAverage
{
public:
    EnsembleAverage (const S& s_)
    : s(s_) 
    {}
    
    double operator() (double beta, const typename S::SMatrix& O) const
    {
        double sum = 0;
        for (size_t i=0; i<s.eigenvalues.size(); i++)
            sum += 
                exp(-beta * s.eigenvalues(i)) * 
                s.eigenvectors.cols(i).transpose() * O * s.eigenvectors.cols(i);
        return sum / partitionFunction(beta);
    }

    //TODO: pull out first/smallest eigenvalue
    double partitionFunction (double beta)
    {
        double Z = 0;
        for (size_t i=0; i<s.eigenvalues.size(); i++)
            Z += std::exp(-beta * s.eigenvalues(i));
        return Z;
    }

private:
    const S& s;
};

#endif // __DIAGONALISER_HPP__ 
