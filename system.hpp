#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>
#include "basis.hpp"
#include <Eigen/Sparse>
#include <Eigen/Dense>
using namespace Eigen;


typedef SparseMatrix<double> SMatrix;
typedef MatrixXd DMatrix;
typedef VectorXd DVector;

enum Spin {UP=0, DOWN=1};

namespace Filter
{
    class SelectAll
    {
    public:
        bool operator() (const State&) const {return true;}
    };
};

namespace Sorter
{
    class DontSort
    {
    public:
        bool operator() (const State&, const State&) const {return false;}
    };
};

template<class System>
class Creator
{
public:
    Creator (const System& s_)
    : s(s_), cs(s.N_orbital)
    {
        for (size_t i=0; i<s.N_orbital; i++)
            constructMatrix(i);
    }

    const SMatrix& operator() (size_t i) const
    {
        return cs[i];
    }

private:
    void constructMatrix (size_t i)
    {
        cs[i] = SMatrix(s.basis.size(), s.basis.size());
        // we expect one entry per column
        cs[i].reserve(s.basis.size());
        for (size_t col=0; col<s.basis.size(); col++)
        {
            cs[i].startVec(col);
            if (s.basis(col)[i] == 1)
                continue;
            State state(s.basis(col));
            state[i] = 1;
            const size_t row = s.basis(state);
            const size_t sum = state.count(i);
            const double sign = (sum%2==0)?1:-1; // probably faster than using (-1)^sum
            cs[i].insertBack(row, col) = sign;
        }
        cs[i].finalize();
    }

    const System& s;
    std::vector<SMatrix> cs;
};

template<class System>
class Annihilator
{
public:
    Annihilator (const System& s_)
    : s(s_), as(s.N_orbital)
    {
        for (size_t i=0; i<s.N_orbital; i++)
            as[i] = SMatrix(s.creator(i).transpose());
    }

    const SMatrix& operator() (size_t i) const
    {
        return as[i];
    }

private:
    const System& s;
    std::vector<SMatrix> as;
};

template<class System>
class Hamiltonian
{
public:
    Hamiltonian (const System& s_)
    : s(s_), H(s.basis.size(), s.basis.size())
    {}

    /**
     * Construct the Hamiltonian matrix.
     *
     * Has to be overwritten by deriving classes.
     */
    void construct() {}

    void diagonalise ()
    {
        /*
         * TODO
         * This assertion is somewhat arbitrary, just because my computers tend
         * to have 4GB ~ 4 10^9 bytes of memory.
         * Also, this is not the solution for the problem. I don't understand
         * why Eigen gives a segfault instead of an exception.
         */
        assert(4E9 > s.basis.size()*s.basis.size()*sizeof(double));
        DMatrix m(H);
        SelfAdjointEigenSolver<DMatrix> es(m);
        eigenvalues = es.eigenvalues();
        eigenvectors = es.eigenvectors();
        Emin = eigenvalues.minCoeff();
    }

    DMatrix eigenvectors;
    DVector eigenvalues;
    double Emin;

protected:
    const System& s;
    SMatrix H;
};

template<class System>
class EnsembleAverage
{
public:
    EnsembleAverage (const System& s_)
    : s(s_)
    {}

    double operator() (double beta, const SMatrix& O) const
    {
        double sum = 0;
        for (int i=0; i<s.H.eigenvalues.size(); i++)
            sum += 
                std::exp(-beta * (s.H.eigenvalues(i) - s.H.Emin)) * 
                s.H.eigenvectors.col(i).adjoint() * O * s.H.eigenvectors.col(i);
        return sum / partitionFunction(beta);
    }

    double partitionFunction (double beta) const
    {
        double Z = 0;
        for (int i=0; i<s.H.eigenvalues.size(); i++)
            Z += std::exp(-beta * (s.H.eigenvalues(i) - s.H.Emin));
        return Z;
    }

private:
    const System& s;
};

template<class Filter, class Sorter>
class MinimalSystem
{
public:
    MinimalSystem (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : N_orbital(N_orbital_), basis(N_orbital, filter, sorter)
    {}

    size_t N_orbital;
    Basis<Filter, Sorter> basis;
};

template<class Filter, class Sorter>
class BasicSystem : public MinimalSystem<Filter, Sorter>
{
public:
    BasicSystem (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : MinimalSystem<Filter, Sorter>(N_orbital_, filter, sorter), creator(*this), 
      annihilator(*this)
    {}

    Creator<BasicSystem> creator;
    Annihilator<BasicSystem> annihilator;
};

/*
template<class Filter, class Sorter>
class System
{
public:
    System (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : N_orbital(N_orbital_), basis(N_orbital, filter, sorter), creator(*this), 
      annihilator(*this)
    {}

    size_t N_orbital;
    Basis<Filter, Sorter> basis;
    Creator<System> creator;
    Annihilator<System> annihilator;
};
*/

#endif // __TEST_HPP__
