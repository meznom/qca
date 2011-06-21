/**
 * @file system.hpp
 * Exact diagonalizer base classes.
 * @author Burkhard Ritter
 * @version 
 * @date 2010-12-10
 */

#ifndef __TEST_HPP__
#define __TEST_HPP__

#include <iostream>
#include "eigenHelpers.hpp"
#include "basis.hpp"

enum Spin {UP=0, DOWN=1};

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
 * H.diagonalize(); //here the real work happens
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
     * Diagonalizes the Hamiltonian matrix.
     *
     * Uses a dense matrix. Relies on Eigen for the eigenvalue decomposition.
     */
    void diagonalize ()
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
        
        eigenvectors = EigenHelpers::denseToSparseBlock(denseEigenvectors, 0, 0, H.rows(), H.rows());
    }

    void diagonalizeBlockWise ()
    {
        eigenvalues = DVector(H.rows());
        eigenvalues.setZero();
        eigenvectors = SMatrix(H.rows(), H.rows());
        eigenvectors.setZero();
        const std::vector<Range>& rs = s.basis.getRanges();
        //we rely on the fact that the ranges are strictly ordered, i.e. a of
        //the current range is >= b of the previous range
        for (size_t i=0; i<rs.size(); i++)
        {
            const int a = rs[i].a;
            const int b = rs[i].b;
            DMatrix m = EigenHelpers::sparseToDenseBlock(H, a, a, b-a, b-a);
            SelfAdjointEigenSolver<DMatrix> es(m);
            DVector blockEigenvalues = es.eigenvalues();
            DMatrix blockEigenvectors = es.eigenvectors();
            eigenvalues.segment(a, b-a) = blockEigenvalues;
            denseToSparseBlockInPlace(blockEigenvectors, eigenvectors, a, a, b-a, b-a);
        }
        eigenvectors.finalize();
        Emin = eigenvalues.minCoeff();
    }

    SMatrix eigenvectors;
    DVector eigenvalues;
    double Emin;

protected:
    void denseToSparseBlockInPlace (const DMatrix& dm, SMatrix& sm, int i, int j, int p, int q)
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
        //for (int i=0; i<s.H.eigenvalues.size(); i++)
        //    sum += 
        //        std::exp(-beta * (s.H.eigenvalues(i) - s.H.Emin)) * 
        //        s.H.eigenvectors.col(i).adjoint() * O * s.H.eigenvectors.col(i);
       
        //TODO: check how efficient this is
        for (int i=0; i<s.H.eigenvalues.size(); i++)
        {
            SMatrix m = s.H.eigenvectors.col(i).adjoint() * O * s.H.eigenvectors.col(i);
            assert(m.size() == 1);
            if (m.nonZeros() != 0)
                sum += std::exp(-beta * (s.H.eigenvalues(i) - s.H.Emin)) * m.coeffRef(0,0);
        }
        
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
 * Only contains the basis. 
 *
 * Typical usage is to subclass and specify a custom Hamiltonian (subclassed
 * from Hamiltonian) as mySystem.H and probably some observables of interest
 * (say the particle number) as well as the EnsembleAverage.
 */
class MinimalSystem
{
public:
    /**
     * Construct the minimal system. 
     *
     * @param N_orbital number of orbitals
     */
    MinimalSystem (size_t N_orbital_)
    : N_orbital(N_orbital_), basis(N_orbital)
    {}

    void construct ()
    {
        basis.construct();
    }

    size_t N_orbital;
    Basis basis;
};

/**
 * Basic base class for fermionic quantum systems.
 *
 * In addition to the basis the basic system conveniently defines creator and
 * annihilator.
 */
class BasicSystem : public MinimalSystem
{
public:
    /**
     * Construct the basic system.
     *
     * @param N_orbital_
     */
    BasicSystem (size_t N_orbital_)
    : MinimalSystem(N_orbital_), creator(*this), annihilator(*this)
    {}

    void construct ()
    {
        MinimalSystem::construct();
        creator.construct();
        annihilator.construct();
    }

    Creator<BasicSystem> creator;
    Annihilator<BasicSystem> annihilator;
};

#endif // __TEST_HPP__
