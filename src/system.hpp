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
public:
    Hamiltonian (System& s_)
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
     *
     * Eigenvectors are stored in a sparse matrix so the eigenvalues() and
     * eigenvectors() accessors should be used for best performance.
     */
    void diagonalizeNoSymmetries ()
    {
        dEigenvalues.resize(0);
        dEigenvectors.resize(0);

        /*
         * This assertion is somewhat arbitrary, just because my computers tend
         * to have 4GB ~ 4 10^9 bytes of memory.
         * Also, this is not the solution for the problem. I don't understand
         * why Eigen gives a segfault instead of an exception.
         */
        assert(4E9 > s.basis.size()*s.basis.size()*sizeof(double));
        SelfAdjointEigenSolver<DMatrix> es(H);
        sEigenvalues = es.eigenvalues();
        minEnergy = sEigenvalues.minCoeff();
        const DMatrix& denseEigenvectors = es.eigenvectors();
        sEigenvectors = denseToSparseBlock(denseEigenvectors, 0, 0, H.rows(), H.rows());
    }

    /**
     * Diagonalizes the Hamiltonian matrix using symmetries.
     *
     * Eploits symmetries defined in the basis (i.e. SymmetryOperators). The
     * symmetries cast the Hamiltonian in a block-diagonal form and this method
     * then diagonalizes blockwise. Usually the method is faster than the
     * simpler diagonalizeNoSymmetries.
     *
     * Eigenvectors are stored in a sparse matrix so the eigenvalues() and
     * eigenvectors() accessors should be used for best performance.
     */
    void diagonalizeUsingSymmetries ()
    {
        dEigenvalues.resize(0);
        dEigenvectors.resize(0);

        sEigenvalues = DVector(H.rows());
        sEigenvalues.setZero();
        sEigenvectors = SMatrix(H.rows(), H.rows());
        sEigenvectors.setZero();

        SelfAdjointEigenSolver<DMatrix> es;

        /* 
         * sorting the ranges according to size means we reduce the number of
         * dynamic memory allocations in the EigenSolver. In practice the
         * performance improvement is negligible. Starting with largest blocks
         * first also means we could improve overall memory usage (by allocating
         * memory for the eigenvectors only as they grow).
         */
        std::vector<Range> rs = s.basis.getRanges();
        std::sort(rs.begin(), rs.end(), CompareSizeOfRanges());
        size_t evsSize = 0;
        for (size_t i=0; i<rs.size(); i++)
            evsSize += (rs[i].b-rs[i].a)*(rs[i].b-rs[i].a);
        sEigenvectors.reserve(evsSize);
            
        size_t index=0;
        for (size_t i=0; i<rs.size(); i++)
        {
            const int a = rs[i].a;
            const int b = rs[i].b;
            //useful debug output when diagonalizing very large Hamiltonians
            //std::cerr << "-> " << "Diagonalizing range " << i << " out of " 
            //          << rs.size() << " ranges." << std::endl
            //          << "-> Size of range " << i << ": " << b-a << std::endl;
            assert(4E9 > (b-a)*(b-a)*sizeof(double));
            es.compute(H.block(a,a,b-a,b-a));
            sEigenvalues.segment(index, b-a) = es.eigenvalues();
            const DMatrix& blockEigenvectors = es.eigenvectors();
            denseToSparseBlockInPlace(blockEigenvectors, sEigenvectors, a, index, b-a, b-a);
            index += b-a;
        }
        sEigenvectors.finalize();
        minEnergy = sEigenvalues.minCoeff();
    }

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
    void diagonalizeUsingSymmetriesBySectors ()
    {
        sEigenvalues.resize(0);
        sEigenvectors.resize(0,0);
        const std::vector<Range>& rs = s.basis.getRanges();
        dEigenvalues.resize(rs.size());
        dEigenvectors.resize(rs.size());
        std::vector<size_t> is(rs.size());
        for (size_t i=0; i<is.size(); i++)
            is[i] = i;
        std::sort(is.begin(), is.end(), SortIndicesAccordingToSizeOfRanges(rs));

        const size_t n_largest = rs[is[0]].b-rs[is[0]].a;
        SelfAdjointEigenSolver<DMatrix> es(n_largest*n_largest);
        
        for (size_t i=0; i<is.size(); i++)
        {
            const int a = rs[is[i]].a;
            const int b = rs[is[i]].b;
            // useful debug output when diagonalizing very large Hamiltonians
            // std::cerr << "-> " << "Diagonalizing range " << is[i] << " out of " 
            //           << rs.size() << " ranges." << std::endl
            //           << "-> Size of range " << is[i] << ": " << b-a << std::endl;
            assert(4E9 > (b-a)*(b-a)*sizeof(double));
            es.compute(H.block(a,a,b-a,b-a));
            dEigenvalues[is[i]] = es.eigenvalues();
            dEigenvectors[is[i]] = es.eigenvectors();
        }
        minEnergy = dEigenvalues[0].minCoeff();
        for (size_t i=1; i<dEigenvalues.size(); i++)
            if (dEigenvalues[i].minCoeff() < minEnergy) 
                minEnergy = dEigenvalues[i].minCoeff();
    }

    void diagonalize ()
    {
        diagonalizeUsingSymmetriesBySectors();
    }

    const DVector& eigenvalues ()
    {
        if (sEigenvalues.size() == 0 && dEigenvalues.size() > 0)
        {
            sEigenvalues.resize(H.rows());
            int k=0;
            for (size_t i=0; i<dEigenvalues.size(); i++)
                for (int j=0; j<dEigenvalues[i].size(); j++)
                {
                    assert(k<sEigenvalues.size());
                    sEigenvalues(k) = dEigenvalues[i](j);
                    k++;
                }
        }
        return sEigenvalues;
    }

    const SMatrix& eigenvectors () 
    {
        if (sEigenvectors.size() == 0 && dEigenvectors.size() > 0)
        {
            sEigenvectors = SMatrix(H.rows(), H.rows());
            sEigenvectors.setZero();
            int k=0;
            for (size_t i=0; i<dEigenvectors.size(); i++)
            {
                const int size = dEigenvectors[i].rows();
                denseToSparseBlockInPlace(dEigenvectors[i], sEigenvectors, k, k, size, size);
                k += size;
            }
            sEigenvectors.finalize();
        }
        return sEigenvectors;
    }

    const std::vector<DVector>& eigenvaluesBySectors ()
    {
        if (dEigenvalues.size() == 0 && sEigenvalues.size() > 0)
            dEigenvalues.push_back(sEigenvalues);
        return dEigenvalues;
    }

    const std::vector<DMatrix>& eigenvectorsBySectors ()
    {
        if (dEigenvectors.size() == 0 && sEigenvectors.size() > 0)
            dEigenvectors.push_back(sEigenvectors);
        return dEigenvectors;
    }

    double Emin () const
    {
        return minEnergy;
    }

protected:
    SMatrix denseToSparseBlock (const DMatrix& dm, int i, int j, int p, int q)
    {
        SMatrix sm(p,q);
        for (int l=0; l<q; l++)
        {
            sm.startVec(l);
            for (int k=0; k<p; k++)
                if (dm(k+i,l+j) != 0)
                    sm.insertBack(k,l) = dm(k+i,l+j);
        }
        sm.finalize();
        return sm;
    }

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

    class CompareSizeOfRanges
    {
    public:
        bool operator() (const Range& r1, const Range& r2) const
        {
            return (r1.b-r1.a)>(r2.b-r2.a);
        }
    };

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
    SMatrix sEigenvectors; //simple or sparse eigenvectors
    DVector sEigenvalues; //simple eigenvalues
    std::vector<DMatrix> dEigenvectors; //dense eigenvectors
    std::vector<DVector> dEigenvalues; //dense eigenvalues
    double minEnergy;
};

/**
 * Calculate the ensemble average.
 *
 * Accesses eigenvectors as one big sparse matrix. Thus it should be used
 * whenever Hamiltonian::diagonalizeNoSymmetries or
 * Hamiltonian::diagonalizeUsingSymmetries is used, for best performance.
 * Example usage:
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
        const DVector& eigenvalues = s.H.eigenvalues();
        const SMatrix& eigenvectors = s.H.eigenvectors();
        for (int i=0; i<eigenvalues.size(); i++)
        {
            const SMatrix& bra = eigenvectors.col(i).adjoint();
            const SMatrix& ket = eigenvectors.col(i);
            SMatrix m = bra * O * ket;
            // Does not work. Should be a bug in Eigen.
            //SMatrix m = eigenvectors.col(i).adjoint() * O * eigenvectors.col(i);
            assert(m.size() == 1);
            if (m.nonZeros() != 0)
                sum += std::exp(-beta * (eigenvalues(i) - s.H.Emin())) * m.coeffRef(0,0);
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
        const DVector& eigenvalues = s.H.eigenvalues();
        for (int i=0; i<eigenvalues.size(); i++)
            Z += std::exp(-beta * (eigenvalues(i) - s.H.Emin()));
        return Z;
    }

private:
    System& s;
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
class EnsembleAverageBySectors
{
public:
    EnsembleAverageBySectors (System& s_)
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

#endif // __TEST_HPP__
