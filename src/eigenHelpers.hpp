/**
 * @file eigenHelpers.hpp
 * Helper and supplementary functions for Eigen.
 * @author Burkhard Ritter
 * @version 
 * @date 2011-06-20
 */

#ifndef __EIGENHELPERS_HPP__
#define __EIGENHELPERS_HPP__

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
using namespace Eigen;

typedef SparseMatrix<double> SMatrix;
typedef MatrixXd DMatrix;
typedef VectorXd DVector;

namespace EigenHelpers
{

// TODO: update to use Eigen's new SparseMatrix interface (not startVec, insertBack
// anymore)

/**
 * Block function for sparse matrices, returning a dense matrix. 
 * 
 * Eigen doesn't implement block() for sparse matrices, only for dense
 * matrices. Hence this small helper function. Can also be used to convert
 * sparse to dense matrices.
 */
inline DMatrix sparseToDenseBlock (const SMatrix& sm, int i, int j, int p, int q)
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

inline void sparseToDenseBlockInPlace (const SMatrix& sm, DMatrix& dm, int i, int j, int p, int q)
{
    dm.resize(p,q);
    dm.setZero();
    for (int l=0; l<q; l++)
    {
        for (SMatrix::InnerIterator it(sm,l+j); it; ++it)
            if (it.row()>=i && it.row()<i+p)
                dm(it.row()-i,it.col()-j) = it.value();
    }
}

/**
 * Block function for dense matrices, returning a sparse matrix. 
 * 
 * Can also be used to convert dense to sparse matrices.
 */
inline SMatrix denseToSparseBlock (const DMatrix& dm, int i, int j, int p, int q)
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

/**
 * Block function for sparse matrices, returning a sparse matrix. 
 * 
 * Eigen doesn't implement block() for sparse matrices, only for dense
 * matrices. Hence this small helper function.
 */
inline SMatrix sparseToSparseBlock (const SMatrix& om, int i, int j, int p, int q)
{
    SMatrix nm(p,q);
    for (int l=0; l<q; l++)
    {
        nm.startVec(l);
        for (SMatrix::InnerIterator it(om,l+j); it; ++it)
            if (it.row()>=i && it.row()<i+p)
                nm.insertBack(it.row()-i,it.col()-j) = it.value();
    }
    nm.finalize();
    return nm;
}

template<typename _MatrixType>
class MySelfAdjointEigenSolver : public SelfAdjointEigenSolver<_MatrixType>
{
private:
    typedef MySelfAdjointEigenSolver<_MatrixType> Self;
    typedef SelfAdjointEigenSolver<_MatrixType> Base;
public:
    MySelfAdjointEigenSolver()
        : Base()
    {}

    MySelfAdjointEigenSolver(const _MatrixType& matrix, int options = ComputeEigenvectors)
        : Base(matrix, options)
    {}

    MySelfAdjointEigenSolver& computeBlock(const SMatrix& matrix, 
                                           int i, int j, int p, int q, 
                                           int options = ComputeEigenvectors);
};

/*
 * This is copied from Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h
 * SelfAdjointEigenSolver::compute and only slightly modified to allow the
 * diagonalization of a block.
 */
template<typename MatrixType>
MySelfAdjointEigenSolver<MatrixType>& MySelfAdjointEigenSolver<MatrixType>
::computeBlock(const SMatrix& matrix, int i, int j, int p, int q, int options)
{
  eigen_assert(p == q);
  eigen_assert((options&~(EigVecMask|GenEigMask))==0
          && (options&EigVecMask)!=EigVecMask
          && "invalid option parameter");
  bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
  typename Self::Index n = p;
  Self::m_eivalues.resize(n,1);

  if(n==1)
  {
    Self::m_eivalues.coeffRef(0,0) = internal::real(matrix.coeff(i,j));
    if(computeEigenvectors)
      Self::m_eivec.setOnes(n,n);
    Self::m_info = Success;
    Self::m_isInitialized = true;
    Self::m_eigenvectorsOk = computeEigenvectors;
    return *this;
  }

  // declare some aliases
  typename Self::RealVectorType& diag = Self::m_eivalues;
  MatrixType& mat = Self::m_eivec;

  // map the matrix coefficients to [-1:1] to avoid over- and underflow.
  sparseToDenseBlockInPlace(matrix, mat, i, j, p, q);
  typename Self::RealScalar scale = mat.cwiseAbs().maxCoeff();
  if(scale==typename Self::RealScalar(0)) scale = typename Self::RealScalar(1);
  mat /= scale;
  Self::m_subdiag.resize(n-1);
  internal::tridiagonalization_inplace(mat, diag, Self::m_subdiag, computeEigenvectors);
  
  typename Self::Index end = n-1;
  typename Self::Index start = 0;
  typename Self::Index iter = 0; // total number of iterations

  while (end>0)
  {
    for (typename Self::Index i = start; i<end; ++i)
      if (internal::isMuchSmallerThan(internal::abs(Self::m_subdiag[i]),(internal::abs(diag[i])+internal::abs(diag[i+1]))))
        Self::m_subdiag[i] = 0;

    // find the largest unreduced block
    while (end>0 && Self::m_subdiag[end-1]==0)
    {
      end--;
    }
    if (end<=0)
      break;

    // if we spent too many iterations, we give up
    iter++;
    if(iter > Self::m_maxIterations * n) break;

    start = end - 1;
    while (start>0 && Self::m_subdiag[start-1]!=0)
      start--;

    internal::tridiagonal_qr_step<MatrixType::Flags&RowMajorBit ? RowMajor : ColMajor>(diag.data(), Self::m_subdiag.data(), start, end, computeEigenvectors ? Self::m_eivec.data() : (typename Self::Scalar*)0, n);
  }

  if (iter <= Self::m_maxIterations * n)
    Self::m_info = Success;
  else
    Self::m_info = NoConvergence;

  // Sort eigenvalues and corresponding vectors.
  // TODO make the sort optional ?
  // TODO use a better sort algorithm !!
  if (Self::m_info == Success)
  {
    for (typename Self::Index i = 0; i < n-1; ++i)
    {
      typename Self::Index k;
      Self::m_eivalues.segment(i,n-i).minCoeff(&k);
      if (k > 0)
      {
        std::swap(Self::m_eivalues[i], Self::m_eivalues[k+i]);
        if(computeEigenvectors)
          Self::m_eivec.col(i).swap(Self::m_eivec.col(k+i));
      }
    }
  }
  
  // scale back the eigen values
  Self::m_eivalues *= scale;

  Self::m_isInitialized = true;
  Self::m_eigenvectorsOk = computeEigenvectors;
  return *this;
}

};

#endif // __EIGENHELPERS_HPP__
