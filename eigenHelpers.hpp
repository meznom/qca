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
using namespace Eigen;

typedef SparseMatrix<double> SMatrix;
typedef MatrixXd DMatrix;
typedef VectorXd DVector;

namespace EigenHelpers
{

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


/*
 * This is copied from Eigen/src/Eigenvalues/SelfAdjointEigenSolver.h and only
 * slightly modified.
 */


#include "Eigen/src/Eigenvalues/EigenvaluesCommon.h"
#include "Eigen/src/Eigenvalues/Tridiagonalization.h"

/** \eigenvalues_module \ingroup Eigenvalues_Module
  *
  *
  * \class MySelfAdjointEigenSolver
  *
  * \brief Computes eigenvalues and eigenvectors of selfadjoint matrices
  *
  * \tparam _MatrixType the type of the matrix of which we are computing the
  * eigendecomposition; this is expected to be an instantiation of the Matrix
  * class template.
  *
  * A matrix \f$ A \f$ is selfadjoint if it equals its adjoint. For real
  * matrices, this means that the matrix is symmetric: it equals its
  * transpose. This class computes the eigenvalues and eigenvectors of a
  * selfadjoint matrix. These are the scalars \f$ \lambda \f$ and vectors
  * \f$ v \f$ such that \f$ Av = \lambda v \f$.  The eigenvalues of a
  * selfadjoint matrix are always real. If \f$ D \f$ is a diagonal matrix with
  * the eigenvalues on the diagonal, and \f$ V \f$ is a matrix with the
  * eigenvectors as its columns, then \f$ A = V D V^{-1} \f$ (for selfadjoint
  * matrices, the matrix \f$ V \f$ is always invertible). This is called the
  * eigendecomposition.
  *
  * The algorithm exploits the fact that the matrix is selfadjoint, making it
  * faster and more accurate than the general purpose eigenvalue algorithms
  * implemented in EigenSolver and ComplexEigenSolver.
  *
  * Only the \b lower \b triangular \b part of the input matrix is referenced.
  *
  * Call the function compute() to compute the eigenvalues and eigenvectors of
  * a given matrix. Alternatively, you can use the
  * MySelfAdjointEigenSolver(const MatrixType&, int) constructor which computes
  * the eigenvalues and eigenvectors at construction time. Once the eigenvalue
  * and eigenvectors are computed, they can be retrieved with the eigenvalues()
  * and eigenvectors() functions.
  *
  * The documentation for MySelfAdjointEigenSolver(const MatrixType&, int)
  * contains an example of the typical use of this class.
  *
  * To solve the \em generalized eigenvalue problem \f$ Av = \lambda Bv \f$ and
  * the likes, see the class GeneralizedSelfAdjointEigenSolver.
  *
  * \sa MatrixBase::eigenvalues(), class EigenSolver, class ComplexEigenSolver
  */
template<typename _MatrixType> class MySelfAdjointEigenSolver
{
  public:

    typedef _MatrixType MatrixType;
    enum {
      Size = MatrixType::RowsAtCompileTime,
      ColsAtCompileTime = MatrixType::ColsAtCompileTime,
      Options = MatrixType::Options,
      MaxColsAtCompileTime = MatrixType::MaxColsAtCompileTime
    };

    /** \brief Scalar type for matrices of type \p _MatrixType. */
    typedef typename MatrixType::Scalar Scalar;
    typedef typename MatrixType::Index Index;

    /** \brief Real scalar type for \p _MatrixType.
      *
      * This is just \c Scalar if #Scalar is real (e.g., \c float or
      * \c double), and the type of the real part of \c Scalar if #Scalar is
      * complex.
      */
    typedef typename NumTraits<Scalar>::Real RealScalar;

    /** \brief Type for vector of eigenvalues as returned by eigenvalues().
      *
      * This is a column vector with entries of type #RealScalar.
      * The length of the vector is the size of \p _MatrixType.
      */
    typedef typename internal::plain_col_type<MatrixType, RealScalar>::type RealVectorType;
    typedef Tridiagonalization<MatrixType> TridiagonalizationType;

    /** \brief Default constructor for fixed-size matrices.
      *
      * The default constructor is useful in cases in which the user intends to
      * perform decompositions via compute(). This constructor
      * can only be used if \p _MatrixType is a fixed-size matrix; use
      * MySelfAdjointEigenSolver(Index) for dynamic-size matrices.
      *
      * Example: \include SelfAdjointEigenSolver_SelfAdjointEigenSolver.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_SelfAdjointEigenSolver.out
      */
    MySelfAdjointEigenSolver()
        : m_eivec(),
          m_eivalues(),
          m_subdiag(),
          m_isInitialized(false)
    { }

    /** \brief Constructor, pre-allocates memory for dynamic-size matrices.
      *
      * \param [in]  size  Positive integer, size of the matrix whose
      * eigenvalues and eigenvectors will be computed.
      *
      * This constructor is useful for dynamic-size matrices, when the user
      * intends to perform decompositions via compute(). The \p size
      * parameter is only used as a hint. It is not an error to give a wrong
      * \p size, but it may impair performance.
      *
      * \sa compute() for an example
      */
    MySelfAdjointEigenSolver(Index size)
        : m_eivec(size, size),
          m_eivalues(size),
          m_subdiag(size > 1 ? size - 1 : 1),
          m_isInitialized(false)
    {}

    /** \brief Constructor; computes eigendecomposition of given matrix.
      *
      * \param[in]  matrix  Selfadjoint matrix whose eigendecomposition is to
      *    be computed. Only the lower triangular part of the matrix is referenced.
      * \param[in]  options Can be ComputeEigenvectors (default) or EigenvaluesOnly.
      *
      * This constructor calls compute(const MatrixType&, int) to compute the
      * eigenvalues of the matrix \p matrix. The eigenvectors are computed if
      * \p options equals ComputeEigenvectors.
      *
      * Example: \include SelfAdjointEigenSolver_SelfAdjointEigenSolver_MatrixType.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_SelfAdjointEigenSolver_MatrixType.out
      *
      * \sa compute(const MatrixType&, int)
      */
    MySelfAdjointEigenSolver(const SMatrix& matrix, int options = ComputeEigenvectors)
      : m_eivec(matrix.rows(), matrix.cols()),
        m_eivalues(matrix.cols()),
        m_subdiag(matrix.rows() > 1 ? matrix.rows() - 1 : 1),
        m_isInitialized(false)
    {
      compute(matrix, options);
    }

    MySelfAdjointEigenSolver& computeBlock(const SMatrix& matrix, 
                                           int i, int j, int p, int q, 
                                           int options = ComputeEigenvectors);

    /** \brief Computes eigendecomposition of given matrix.
      *
      * \param[in]  matrix  Selfadjoint matrix whose eigendecomposition is to
      *    be computed. Only the lower triangular part of the matrix is referenced.
      * \param[in]  options Can be ComputeEigenvectors (default) or EigenvaluesOnly.
      * \returns    Reference to \c *this
      *
      * This function computes the eigenvalues of \p matrix.  The eigenvalues()
      * function can be used to retrieve them.  If \p options equals ComputeEigenvectors,
      * then the eigenvectors are also computed and can be retrieved by
      * calling eigenvectors().
      *
      * This implementation uses a symmetric QR algorithm. The matrix is first
      * reduced to tridiagonal form using the Tridiagonalization class. The
      * tridiagonal matrix is then brought to diagonal form with implicit
      * symmetric QR steps with Wilkinson shift. Details can be found in
      * Section 8.3 of Golub \& Van Loan, <i>%Matrix Computations</i>.
      *
      * The cost of the computation is about \f$ 9n^3 \f$ if the eigenvectors
      * are required and \f$ 4n^3/3 \f$ if they are not required.
      *
      * This method reuses the memory in the MySelfAdjointEigenSolver object that
      * was allocated when the object was constructed, if the size of the
      * matrix does not change.
      *
      * Example: \include SelfAdjointEigenSolver_compute_MatrixType.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_compute_MatrixType.out
      *
      * \sa SelfAdjointEigenSolver(const MatrixType&, int)
      */
    //MySelfAdjointEigenSolver& compute(const MatrixType& matrix, int options = ComputeEigenvectors);
    //template<typename Derived>
    //MySelfAdjointEigenSolver& compute(const MatrixBase<Derived>& matrix, int options = ComputeEigenvectors);
    MySelfAdjointEigenSolver& compute(const SMatrix& matrix, int options = ComputeEigenvectors)
    {
        return computeBlock(matrix, 0, 0, matrix.rows(), matrix.cols(), options);
    }

    /** \brief Returns the eigenvectors of given matrix.
      *
      * \returns  A const reference to the matrix whose columns are the eigenvectors.
      *
      * \pre The eigenvectors have been computed before.
      *
      * Column \f$ k \f$ of the returned matrix is an eigenvector corresponding
      * to eigenvalue number \f$ k \f$ as returned by eigenvalues().  The
      * eigenvectors are normalized to have (Euclidean) norm equal to one. If
      * this object was used to solve the eigenproblem for the selfadjoint
      * matrix \f$ A \f$, then the matrix returned by this function is the
      * matrix \f$ V \f$ in the eigendecomposition \f$ A = V D V^{-1} \f$.
      *
      * Example: \include SelfAdjointEigenSolver_eigenvectors.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_eigenvectors.out
      *
      * \sa eigenvalues()
      */
    const MatrixType& eigenvectors() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec;
    }

    /** \brief Returns the eigenvalues of given matrix.
      *
      * \returns A const reference to the column vector containing the eigenvalues.
      *
      * \pre The eigenvalues have been computed before.
      *
      * The eigenvalues are repeated according to their algebraic multiplicity,
      * so there are as many eigenvalues as rows in the matrix. The eigenvalues
      * are sorted in increasing order.
      *
      * Example: \include SelfAdjointEigenSolver_eigenvalues.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_eigenvalues.out
      *
      * \sa eigenvectors(), MatrixBase::eigenvalues()
      */
    const RealVectorType& eigenvalues() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      return m_eivalues;
    }

    /** \brief Computes the positive-definite square root of the matrix.
      *
      * \returns the positive-definite square root of the matrix
      *
      * \pre The eigenvalues and eigenvectors of a positive-definite matrix
      * have been computed before.
      *
      * The square root of a positive-definite matrix \f$ A \f$ is the
      * positive-definite matrix whose square equals \f$ A \f$. This function
      * uses the eigendecomposition \f$ A = V D V^{-1} \f$ to compute the
      * square root as \f$ A^{1/2} = V D^{1/2} V^{-1} \f$.
      *
      * Example: \include SelfAdjointEigenSolver_operatorSqrt.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_operatorSqrt.out
      *
      * \sa operatorInverseSqrt(),
      *     \ref MatrixFunctions_Module "MatrixFunctions Module"
      */
    MatrixType operatorSqrt() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec * m_eivalues.cwiseSqrt().asDiagonal() * m_eivec.adjoint();
    }

    /** \brief Computes the inverse square root of the matrix.
      *
      * \returns the inverse positive-definite square root of the matrix
      *
      * \pre The eigenvalues and eigenvectors of a positive-definite matrix
      * have been computed before.
      *
      * This function uses the eigendecomposition \f$ A = V D V^{-1} \f$ to
      * compute the inverse square root as \f$ V D^{-1/2} V^{-1} \f$. This is
      * cheaper than first computing the square root with operatorSqrt() and
      * then its inverse with MatrixBase::inverse().
      *
      * Example: \include SelfAdjointEigenSolver_operatorInverseSqrt.cpp
      * Output: \verbinclude SelfAdjointEigenSolver_operatorInverseSqrt.out
      *
      * \sa operatorSqrt(), MatrixBase::inverse(),
      *     \ref MatrixFunctions_Module "MatrixFunctions Module"
      */
    MatrixType operatorInverseSqrt() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      eigen_assert(m_eigenvectorsOk && "The eigenvectors have not been computed together with the eigenvalues.");
      return m_eivec * m_eivalues.cwiseInverse().cwiseSqrt().asDiagonal() * m_eivec.adjoint();
    }

    /** \brief Reports whether previous computation was successful.
      *
      * \returns \c Success if computation was succesful, \c NoConvergence otherwise.
      */
    ComputationInfo info() const
    {
      eigen_assert(m_isInitialized && "SelfAdjointEigenSolver is not initialized.");
      return m_info;
    }

    /** \brief Maximum number of iterations.
      *
      * Maximum number of iterations allowed for an eigenvalue to converge.
      */
    static const int m_maxIterations = 30;

  protected:
    MatrixType m_eivec;
    RealVectorType m_eivalues;
    typename TridiagonalizationType::SubDiagonalType m_subdiag;
    ComputationInfo m_info;
    bool m_isInitialized;
    bool m_eigenvectorsOk;
};

template<typename MatrixType>
//template<typename Derived>
//::compute(const MatrixBase<Derived>& matrix, int options)
MySelfAdjointEigenSolver<MatrixType>& MySelfAdjointEigenSolver<MatrixType>
::computeBlock(const SMatrix& matrix, int i, int j, int p, int q, int options)
{
  eigen_assert(p == q);
  eigen_assert((options&~(EigVecMask|GenEigMask))==0
          && (options&EigVecMask)!=EigVecMask
          && "invalid option parameter");
  bool computeEigenvectors = (options&ComputeEigenvectors)==ComputeEigenvectors;
  Index n = p;
  m_eivalues.resize(n,1);

  if(n==1)
  {
    m_eivalues.coeffRef(0,0) = internal::real(matrix.coeff(i,j));
    if(computeEigenvectors)
    {
      m_eivec.resize(1,1);
      m_eivec.setOnes();
    }
    m_info = Success;
    m_isInitialized = true;
    m_eigenvectorsOk = computeEigenvectors;
    return *this;
  }

  // declare some aliases
  RealVectorType& diag = m_eivalues;
  MatrixType& mat = m_eivec;

  // map the matrix coefficients to [-1:1] to avoid over- and underflow.
  //RealScalar scale = matrix.cwiseAbs().maxCoeff();
  //if(scale==Scalar(0)) scale = 1;
  //mat = matrix / scale;
  //mat = sparseToDenseBlock(matrix, i, j, p, q);
  sparseToDenseBlockInPlace(matrix, mat, i, j, p, q);
  RealScalar scale = mat.cwiseAbs().maxCoeff();
  if(scale==Scalar(0)) scale = 1;
  mat /= scale;

  m_subdiag.resize(n-1);
  internal::tridiagonalization_inplace(mat, diag, m_subdiag, computeEigenvectors);
  
  Index end = n-1;
  Index start = 0;
  Index iter = 0; // number of iterations we are working on one element

  while (end>0)
  {
    for (Index i = start; i<end; ++i)
      if (internal::isMuchSmallerThan(internal::abs(m_subdiag[i]),(internal::abs(diag[i])+internal::abs(diag[i+1]))))
        m_subdiag[i] = 0;

    // find the largest unreduced block
    while (end>0 && m_subdiag[end-1]==0)
    {
      iter = 0;
      end--;
    }
    if (end<=0)
      break;

    // if we spent too many iterations on the current element, we give up
    iter++;
    if(iter > m_maxIterations) break;

    start = end - 1;
    while (start>0 && m_subdiag[start-1]!=0)
      start--;

    internal::tridiagonal_qr_step<MatrixType::Flags&RowMajorBit ? RowMajor : ColMajor>(diag.data(), m_subdiag.data(), start, end, computeEigenvectors ? m_eivec.data() : (Scalar*)0, n);
  }

  if (iter <= m_maxIterations)
    m_info = Success;
  else
    m_info = NoConvergence;

  // Sort eigenvalues and corresponding vectors.
  // TODO make the sort optional ?
  // TODO use a better sort algorithm !!
  if (m_info == Success)
  {
    for (Index i = 0; i < n-1; ++i)
    {
      Index k;
      m_eivalues.segment(i,n-i).minCoeff(&k);
      if (k > 0)
      {
        std::swap(m_eivalues[i], m_eivalues[k+i]);
        if(computeEigenvectors)
          m_eivec.col(i).swap(m_eivec.col(k+i));
      }
    }
  }
  
  // scale back the eigen values
  m_eivalues *= scale;

  m_isInitialized = true;
  m_eigenvectorsOk = computeEigenvectors;
  return *this;
}


};

#endif // __EIGENHELPERS_HPP__
