#ifndef __DIAGONALISER_HPP__
#define __DIAGONALISER_HPP__

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

class System
{
public:
    typedef FermionicState<> State;
    typedef State::num_t num_t;
    typedef SparseMatrix<double> SMatrix;
    typedef MatrixXd DMatrix;
    typedef VectorXd DVector;
    typedef std::vector<State>::const_iterator CBIT;

    System (size_t N_orbital_) : N_orbital(N_orbital_), cs(N_orbital_), as(N_orbital_)
    {
        constructBasis(N_orbital, voidFilter, voidSorter);

        
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

    class FilterSelectAll
    {
    public:
        bool operator() (const State&) {return true;}
    };

    class SorterDontSort
    {
    public:
        bool operator() (const State&, const State&) {return false;}
    };

//private:
    std::vector<State> basis;
    num_t N_basis;
    size_t N_orbital;
    FilterSelectAll voidFilter;
    SorterDontSort voidSorter;
    std::vector<SMatrix> cs, as;
    enum Spin {UP=0, DOWN=1};
};

class AndersonModel : public System
{
public:
    AndersonModel (size_t N_bath_)
        : System(2*N_bath_+ 2), N_bath(N_bath_), H(*this, N_bath_)
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

    class Hamiltonian
    {
    public:
        Hamiltonian (const AndersonModel& s_, size_t N_bath_)
            : s(s_), t(1), mu(0), U(0), V(1), Ed(0), N_bath(N_bath_), Hmatrix(s_.basis.size(), s_.basis.size())
        {}

        SMatrix getMatrix ()
        {
            Hmatrix.setZero();
            for (size_t i = 1; i<N_bath+1; i++)
            {
                Hmatrix += - mu * s.c(i,UP) * s.a(i,UP) - mu * s.c(i,DOWN) * s.a(i,DOWN);
                
                Hmatrix += -t * ( 
                    s.c(i,UP) * s.a((i+1)%N_bath,UP) + s.c((i+1)%N_bath,UP) * s.a(i,UP) +
                    s.c(i,DOWN) * s.a((i+1)%N_bath,DOWN) + s.c((i+1)%N_bath,DOWN) * s.a(i,DOWN)
                );
            }
            Hmatrix += Ed * (s.c(0,UP)*s.a(0,UP) + s.c(0,DOWN)*s.a(0,DOWN));
            Hmatrix += V * (
                s.c(0,UP)*s.a(1,UP) + s.c(1,UP)*s.a(0,UP) + 
                s.c(0,DOWN)*s.a(1,DOWN) + s.c(1,DOWN)*s.a(0,DOWN)
            );

            return Hmatrix;
        }

        void diagonalise ()
        {
            DMatrix m(Hmatrix);
            SelfAdjointEigenSolver<DMatrix> es(m);
            eigenvalues = es.eigenvalues();
            eigenvectors = es.eigenvectors();
        }

        const AndersonModel& s;
        double t, mu, U, V, Ed;
        size_t N_bath;
        SMatrix Hmatrix;
        DMatrix eigenvectors;
        DVector eigenvalues;
    };

    Hamiltonian H;
    size_t N_bath;
};



#endif // __DIAGONALISER_HPP__ 
