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

};

#endif // __EIGENHELPERS_HPP__
