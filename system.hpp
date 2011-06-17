/**
 * @file system.hpp
 * Exact diagonaliser base classes.
 * @author Burkhard Ritter
 * @version 
 * @date 2010-12-10
 */

#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>
#include "basis.hpp"

/**
 * Filters to filter the basis states.
 */
namespace Filter
{
    /**
     * Selects all states.
     *
     * Thus it's a filter that doesn't filter at all. It's the default Filter
     * used.
     */
    class SelectAll
    {
    public:
        bool operator() (const State&) const {return true;}
    };
};

/**
 * Sorters to sort the basis states.
 */
namespace Sorter
{
    /**
     * Does not sort the states.
     *
     * Thus it's a sorter that doesn't sort at all. It's the default Sorter
     * used.
     */
    class DontSort
    {
    public:
        bool operator() (const State&, const State&) const {return false;}
    };
};

/**
 * Caching creator. 
 *
 * Constructs one creator martrix for each orbital of the system and caches
 * them.
 *
 * @tparam System
 */
template<class System>
class Creator
{
public:
    Creator (const System& s_)
    : s(s_), cs(s.N_orbital)
    {}

    void construct ()
    {
        for (size_t i=0; i<s.N_orbital; i++)
            constructMatrix(i);
    }

    /**
     * Creator matrix for the given orbital.
     *
     * \f$c^{\dag}_i\f$
     *
     * @param i orbital of the system
     *
     * @return 
     */
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

/**
 * Caching annihilator.
 *
 * Constructs one annihilator matrix for each orbital of the system and caches
 * them.
 *
 * @tparam System
 */
template<class System>
class Annihilator
{
public:
    Annihilator (const System& s_)
    : s(s_), as(s.N_orbital)
    {}

    void construct ()
    {
        for (size_t i=0; i<s.N_orbital; i++)
            as[i] = SMatrix(s.creator(i).transpose());
    }

    /**
     * Annihilator matrix for the given orbital.
     *
     * \f$c_i\f$
     *
     * @param i orbital of the system
     *
     * @return 
     */
    const SMatrix& operator() (size_t i) const
    {
        return as[i];
    }

private:
    const System& s;
    std::vector<SMatrix> as;
};

/**
 * Hamiltonian of the system.
 *
 * The heart of the fermionic quantum system. Although it is an
 * operator it behaves a little bit different than the other operators, which
 * maybe doesn't entirely come as a surprise. It is
 * not a function object. It has to be subclassed and
 * overwritten for each different physical system. Typical usage:
 * @code
 * Hamiltonian H(mySystem);
 * H.construct();
 * H.diagonalise(); //here the real work happens
 * H.eigenvalues(); //use eigenenergies
 * @endcode
 *
 * @tparam System
 */
template<class System>
class Hamiltonian
{
public:
    Hamiltonian (const System& s_)
    : s(s_), H()
    {}

    /**
     * Construct the Hamiltonian matrix.
     *
     * Has to be overwritten by deriving classes.
     */
    void construct() {}

    /**
     * Diagonalises the Hamiltonian matrix.
     *
     * Uses a dense matrix. Relies on Eigen for the eigenvalue decomposition.
     */
    void diagonalise ()
    {
        /*
         * This assertion is somewhat arbitrary, just because my computers tend
         * to have 4GB ~ 4 10^9 bytes of memory.
         * Also, this is not the solution for the problem. I don't understand
         * why Eigen gives a segfault instead of an exception.
         */
        assert(4E9 > s.basis.size()*s.basis.size()*sizeof(double));
        DMatrix m(H);
        SelfAdjointEigenSolver<DMatrix> es(m);
        eigenvalues = es.eigenvalues();
        Emin = eigenvalues.minCoeff();
        DMatrix denseEigenvectors = es.eigenvectors();
        
        eigenvectors = SMatrix(H.rows(), H.rows());
        eigenvectors.setZero();
        denseBlockToSparseMatrix(denseEigenvectors, eigenvectors, 0, 0, H.rows(), H.rows());
        eigenvectors.finalize();
    }

    void diagonaliseBlockWise ()
    {
        eigenvalues = DVector(H.rows());
        eigenvalues.setZero();
        eigenvectors = SMatrix(H.rows(), H.rows());
        eigenvectors.setZero();
        //TODO: eigenvectors.reserve??
        const std::vector<Range>& rs = s.basis.getRanges();
        int oldA = -1;
        int oldB = -1;
        for (size_t i=0; i<rs.size(); i++)
        {
            const int& a = rs[i].a;
            const int& b = rs[i].b;
            
            assert(oldA < a && oldB < b);
            oldA = a;
            oldB = b;

            DMatrix m = block(H, a, a, b-a, b-a);
            SelfAdjointEigenSolver<DMatrix> es(m);
            DVector blockEigenvalues = es.eigenvalues();
            DMatrix blockEigenvectors = es.eigenvectors();
            eigenvalues.segment(a, b-a) = blockEigenvalues;
            denseBlockToSparseMatrix(blockEigenvectors, eigenvectors, a, a, b-a, b-a);
        }
        eigenvectors.finalize();
        Emin = eigenvalues.minCoeff();
    }

    DMatrix block (const SMatrix& sm, int i, int j, int p, int q) const
    {
        DMatrix dm = DMatrix::Zero(p,q);
        for (int l=0; l<q; l++)
        {
            for (SMatrix::InnerIterator it(sm,l+j); it; ++it)
                if (it.row()>=i && it.row()<i+p)
                    dm(it.row()-i,it.col()-j) = it.value();
        }
        return dm;
    }

    void denseBlockToSparseMatrix (const DMatrix& dm, SMatrix& sm, int i, int j, int p, int q)
    {
        assert(dm.rows() == p);
        assert(dm.cols() == q);

        for (int l=0; l<q; l++)
        {
            sm.startVec(l+j);
            for (int k=0; k<p; k++)
                sm.insertBack(k+i,l+j) = dm(k,l);
        }
    }

    //DMatrix eigenvectors;
    SMatrix eigenvectors;
    DVector eigenvalues;
    double Emin;

protected:
    const System& s;
    SMatrix H;
};

/**
 * Calculate the ensemble average.
 *
 * An operator. Example usage:
 * @code
 * EnsembleAverage ensembleAverage(mySystem);
 * MyFunkyOperator O(mySystem);
 * ensembleAverage(beta, O); //will expect and use mySystem.H
 * @endcode
 *
 * Dependencies: System.basis and System.H (Hamiltonian)
 *
 * @tparam System
 */
template<class System>
class EnsembleAverage
{
public:
    EnsembleAverage (const System& s_)
    : s(s_)
    {}

    /**
     * Calculate the ensemble average for the given operator.
     *
     * @param beta = 1/T (temperature)
     * @param O operator matrix
     *
     * @return 
     */
    double operator() (double beta, const SMatrix& O) const
    {
        double sum = 0;
        for (int i=0; i<s.H.eigenvalues.size(); i++)
            sum += 
                std::exp(-beta * (s.H.eigenvalues(i) - s.H.Emin)) * 
                s.H.eigenvectors.col(i).adjoint() * O * s.H.eigenvectors.col(i);
        return sum / partitionFunction(beta);
    }

    /**
     * Calculate the partition function at the given temperature.
     *
     * @param beta = 1/T (temperature)
     *
     * @return 
     */
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

/**
 * Minimal base class for a fermionic quantum system.
 *
 * Only contains the basis. Filter and Sorter are passed on to the basis
 * and are used to filter and sort the basis states.
 *
 * Typical usage is to subclass and specify a custom Hamiltonian (subclassed
 * from Hamiltonian) as mySystem.H and probably some observables of interest
 * (say the particle number) as well as the EnsembleAverage.
 *
 * @tparam Filter
 * @tparam Sorter
 */
template<class Filter, class Sorter>
class MinimalSystem
{
public:
    /**
     * Construct the minimal system. 
     *
     * @param N_orbital number of orbitals
     * @param filter filter to filter basis states
     * @param sorter sorter to sort basis states
     */
    MinimalSystem (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : N_orbital(N_orbital_), basis(N_orbital, filter, sorter)
    {}

    void construct ()
    {
        basis.construct();
    }

    size_t N_orbital;
    Basis<Filter, Sorter> basis;
};

/**
 * Basic base class for fermionic quantum systems.
 *
 * In addition to the basis the basic system conveniently defines creator and
 * annihilator.
 *
 * @tparam Filter
 * @tparam Sorter
 */
template<class Filter, class Sorter>
class BasicSystem : public MinimalSystem<Filter, Sorter>
{
public:
    /**
     * Construct the basic system.
     *
     * @param N_orbital_
     * @param filter
     * @param sorter
     */
    BasicSystem (size_t N_orbital_, const Filter& filter, const Sorter& sorter)
    : MinimalSystem<Filter, Sorter>(N_orbital_, filter, sorter), creator(*this), 
      annihilator(*this)
    {}

    void construct ()
    {
        MinimalSystem<Filter, Sorter>::construct();
        creator.construct();
        annihilator.construct();
    }

    Creator<BasicSystem> creator;
    Annihilator<BasicSystem> annihilator;
};

#endif // __TEST_HPP__
