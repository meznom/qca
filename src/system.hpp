/**
 * @file system.hpp
 * Exact diagonalizer base classes.
 * @author Burkhard Ritter
 * @version 
 * @date 2010-12-10
 */

#ifndef __SYSTEM_HPP__
#define __SYSTEM_HPP__

#include <iostream>
#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "basis.hpp"
using namespace Eigen;

typedef SparseMatrix<double> SMatrix;
typedef MatrixXd DMatrix;
typedef VectorXd DVector;

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
    : s(s_)
    {}

    void construct ()
    {
        const size_t N_orbital = s.basis.numberOfOrbitals();
        cs = std::vector<SMatrix>(N_orbital);
        for (size_t i=0; i<N_orbital; i++)
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
    const SMatrix& operator() (int i) const
    {
        return cs[i];
    }

private:
    // TODO: update to use Eigen's new SparseMatrix interface
    void constructMatrix (int i)
    {
        cs[i] = SMatrix(s.basis.size(), s.basis.size());
        // we expect one entry per column
        cs[i].reserve(VectorXi::Constant(s.basis.size(),1));
        for (size_t col=0; col<s.basis.size(); col++)
        {
            if (s.basis(col)[i] == 1)
                continue;
            State state(s.basis(col));
            state[i] = 1;
            const size_t row = s.basis(state);
            const size_t sum = state.count(i);
            const double sign = (sum%2==0)?1:-1; // probably faster than using (-1)^sum
            cs[i].insert(row, col) = sign;
        }
        cs[i].makeCompressed();
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
    : s(s_)
    {}

    void construct ()
    {
        const size_t N_orbital = s.basis.numberOfOrbitals();
        as = std::vector<SMatrix>(N_orbital);
        for (size_t i=0; i<N_orbital; i++)
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
    const SMatrix& operator() (int i) const
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
protected:
    class SortIndicesAccordingToSizeOfRanges
    {
    public:
        SortIndicesAccordingToSizeOfRanges (const std::vector<Range>& rs_)
            : rs(rs_) 
        {}

        bool operator() (size_t i, size_t j) const
        {
            return (rs[i].b-rs[i].a)>(rs[j].b-rs[j].a);
        }

    private:
        const std::vector<Range>& rs;
    };

    System& s;
    SMatrix H;
    DVector es; // eigenvalues
    // SMatrix vs; // eigenvalues, unused
    std::vector<DVector> es_s; // eigenvalues, sector-wise
    std::vector<DMatrix> vs_s; // eigenvectors, sector-wise
    double E_min;

public:
    Hamiltonian (System& s_)
    : s(s_), H(), E_min(0)
    {}

    /**
     * Construct the Hamiltonian matrix.
     *
     * Has to be overwritten by deriving classes.
     */
    void construct() {}

    /**
     * Diagonalizes the Hamiltonian matrix using symmetries.
     *
     * Similar to diagonalizeUsingSymmetries, but stores eigenvectors as an
     * std::vector of dense matrices. Consequently, the accessors
     * eigenvaluesBySectors() and eigenvectorsBySectors() should be used for
     * best performance.
     *
     * Uses less memory than diagonalizeUsingSymmetries and is also a little bit
     * faster.
     */
    void diagonalize ()
    {
        es.resize(0);
        const std::vector<Range>& rs = s.basis.getRanges();
        es_s.resize(rs.size());
        vs_s.resize(rs.size());
        std::vector<size_t> is(rs.size());
        for (size_t i=0; i<is.size(); i++)
            is[i] = i;
        std::sort(is.begin(), is.end(), SortIndicesAccordingToSizeOfRanges(rs));

        const size_t n_largest = rs[is[0]].b-rs[is[0]].a;
        SelfAdjointEigenSolver<DMatrix> s(n_largest*n_largest);
        
        for (size_t i=0; i<is.size(); i++)
        {
            const int a = rs[is[i]].a;
            const int b = rs[is[i]].b;
            // useful debug output when diagonalizing very large Hamiltonians
            // std::cerr << "-> " << "Diagonalizing range " << is[i] << " out of " 
            //           << rs.size() << " ranges." << std::endl
            //           << "-> Size of range " << is[i] << ": " << b-a << std::endl;
            assert(4E9 > (b-a)*(b-a)*sizeof(double));
            s.compute(H.block(a,a,b-a,b-a));
            es_s[is[i]] = s.eigenvalues();
            vs_s[is[i]] = s.eigenvectors();
        }
        E_min = es_s[0].minCoeff();
        for (size_t i=1; i<es_s.size(); i++)
            if (es_s[i].minCoeff() < E_min) 
                E_min = es_s[i].minCoeff();
    }

    const std::vector<DVector>& eigenvaluesBySectors ()
    {
        return es_s;
    }

    const std::vector<DMatrix>& eigenvectorsBySectors ()
    {
        return vs_s;
    }

    const DVector& eigenvalues ()
    {
        if (es.size() == 0 && es_s.size() > 0)
        {
            es.resize(H.rows());
            int k=0;
            for (size_t i=0; i<es_s.size(); i++)
                for (int j=0; j<es_s[i].size(); j++)
                {
                    assert(k<es.size());
                    es(k) = es_s[i](j);
                    k++;
                }
        }
        return es;
    }

    double Emin () const
    {
        return E_min;
    }
};

/**
 * Calculate the ensemble average "by sectors".
 *
 * "By sectors" essentially means that the eigenvalues are accessed as an
 * std::vector of dense vectors and the eigenvectors are accessed as an
 * std::vector of dense matrices. Hence this class should be used whenever
 * Hamiltonian::diagonlizeUsingSymmetriesBySectors is used, for best
 * performance. Example usage:
 * @code
 * EnsembleAverageBySectors ensembleAverage(mySystem);
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
    EnsembleAverage (System& s_)
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
        size_t index = 0;
        const std::vector<DVector>& eigenvalues = s.H.eigenvaluesBySectors();
        const std::vector<DMatrix>& eigenvectors = s.H.eigenvectorsBySectors();
        for (size_t i=0; i<eigenvalues.size(); i++)
        {
            // note: usually O is very sparse, so it is essential to use a
            // sparse matrix for the block
            const int size = eigenvalues[i].size();
            const SMatrix& O_block = O.block(index,index,size,size);
            for (int j=0; j<size; j++)
                sum += 
                    std::exp(-beta * (eigenvalues[i](j) - s.H.Emin())) * 
                    eigenvectors[i].col(j).adjoint() * O_block * eigenvectors[i].col(j);
            index += size;
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
        const std::vector<DVector>& eigenvalues = s.H.eigenvaluesBySectors();
        for (size_t i=0; i<eigenvalues.size(); i++)
            for (int j=0; j<eigenvalues[i].size(); j++)
                Z += std::exp(-beta * (eigenvalues[i](j) - s.H.Emin()));
        return Z;
    }

private:
    System& s;
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
    : N_orbital(N_orbital_)
    {}

    size_t N_orbital;
    Basis basis;

protected:
    void construct ()
    {
        basis.construct(N_orbital);
    }
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

    Creator<BasicSystem> creator;
    Annihilator<BasicSystem> annihilator;

protected:
    void construct ()
    {
        MinimalSystem::construct();
        creator.construct();
        annihilator.construct();
    }
};

#endif // __SYSTEM_HPP__
