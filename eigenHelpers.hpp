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
    DMatrix sparseToDenseBlock (const SMatrix& sm, int i, int j, int p, int q)
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

    /**
     * Block function for dense matrices, returning a sparse matrix. 
     * 
     * Can also be used to convert dense to sparse matrices.
     */
    SMatrix denseToSparseBlock (const DMatrix& dm, int i, int j, int p, int q)
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
    SMatrix sparseToSparseBlock (const SMatrix& om, int i, int j, int p, int q)
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
};

#endif // __EIGENHELPERS_HPP__
