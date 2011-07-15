#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE eigenHelpers test
#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>
#include "eigenHelpers.hpp"

using namespace EigenHelpers;

SMatrix constructSparseMatrix ()
{
    DynamicSparseMatrix<double> dsm(5,5);
    int k=1;
    for (int i=0; i<5; i++)
        for (int j=0; j<5; j++)
        {
            dsm.coeffRef(i,j) = k;
            k++;
        }
    SMatrix sm(dsm);
    return sm;
}

bool sparseDenseEqual (SMatrix& sm, const DMatrix& dm)
{
    if (sm.rows() != dm.rows() || sm.cols() != dm.cols())
        return false;
    for (int j=0; j<sm.cols(); j++)
        for (int i=0; i<sm.rows(); i++)
            if (sm.coeffRef(i,j) != dm(i,j))
                return false;
    return true;
}

bool operator== (SMatrix& sm1, SMatrix& sm2)
{
    if (sm1.rows() != sm2.rows() || sm1.cols() != sm2.cols())
        return false;
    for (int j=0; j<sm1.cols(); j++)
        for (int i=0; i<sm1.rows(); i++)
            if (sm1.coeffRef(i,j) != sm2.coeffRef(i,j))
                return false;
    return true;
}

BOOST_AUTO_TEST_CASE ( sparse_to_dense_block )
{
    SMatrix sm = constructSparseMatrix();
    DMatrix dm1(5,5);
    dm1 <<  1, 2, 3, 4, 5,
            6, 7, 8, 9,10,
           11,12,13,14,15,
           16,17,18,19,20,
           21,22,23,24,25;
    BOOST_CHECK (sparseDenseEqual(sm,dm1));
    BOOST_CHECK (sparseToDenseBlock(sm, 0, 0, 5, 5) == dm1);

    DMatrix dm2(3,3);
    dm2 <<  7, 8, 9,
           12,13,14,
           17,18,19;
    BOOST_CHECK (sparseToDenseBlock(sm, 1, 1, 3, 3) == dm2);

    DMatrix dm3(2,3);
    dm3 << 18,19,20,
           23,24,25;
    BOOST_CHECK (sparseToDenseBlock(sm, 3, 2, 2, 3) == dm3);
}

BOOST_AUTO_TEST_CASE ( dense_to_sparse_block )
{
    SMatrix sm1 = constructSparseMatrix();
    DMatrix dm1(5,5);
    dm1 <<  1, 2, 3, 4, 5,
            6, 7, 8, 9,10,
           11,12,13,14,15,
           16,17,18,19,20,
           21,22,23,24,25;
    SMatrix sm2 = denseToSparseBlock(dm1, 0, 0, 5, 5);
    BOOST_CHECK (sm1 == sm2);

    DMatrix dm2(7,7);
    dm2 <<  1, 1, 0, 0, 9, 9, 0,
            1, 1, 2, 3, 4, 5, 0,
            1, 6, 7, 8, 9,10, 0,
            1,11,12,13,14,15, 0,
            1,16,17,18,19,20, 0,
            1,21,22,23,24,25, 0,
            0, 0, 0, 0, 0, 0, 0;
    sm2 = denseToSparseBlock(dm2, 1, 1, 5, 5);
    BOOST_CHECK (sm1 == sm2);

    DMatrix dm3(7,7);
    dm3 <<   1, 0, 0, 9, 9, 0,1,
             0, 0, 0, 0, 0, 0,0,
             1, 2, 3, 4, 5, 0,1,
             6, 7, 8, 9,10, 0,1,
            11,12,13,14,15, 0,1,
            16,17,18,19,20, 0,1,
            21,22,23,24,25, 0,1;
    sm2 = denseToSparseBlock(dm3, 2, 0, 5, 5);
    BOOST_CHECK (sm1 == sm2);

    DMatrix dm4(5,5);
    dm4 <<  0, 2, 0, 4, 0,
            0, 0, 0, 0,10,
           11, 0, 0, 0, 0,
           16, 0, 0,19, 0,
           21, 0, 0, 0, 0;
    SMatrix sm4 = denseToSparseBlock(dm4, 0, 0, 5, 5);
    BOOST_CHECK (sm4.nonZeros() == 7);
}

BOOST_AUTO_TEST_CASE ( sparse_to_sparse_block )
{
    SMatrix sm0 = constructSparseMatrix();
    DMatrix dm1(3,3);
    dm1 <<  1, 2, 3,
            6, 7, 8,
           11,12,13;
    SMatrix sm1 = denseToSparseBlock(dm1, 0, 0, 3, 3);
    SMatrix sm2 = sparseToSparseBlock(sm0, 0, 0, 3, 3);
    BOOST_CHECK (sm1 == sm2);

    DMatrix dm2(3,2);
    dm2 <<  9,10,
           14,15,
           19,20;
    sm1 = denseToSparseBlock(dm2, 0, 0, 3, 2);
    sm2 = sparseToSparseBlock(sm0, 1, 3, 3, 2);
    BOOST_CHECK (sm1 == sm2);
}

BOOST_AUTO_TEST_CASE ( performance_of_MySelfAdjointEigenSolver_computeBlock )
{
#ifdef NDEBUG
    std::clock_t startCPUTime, endCPUTime;
    double cpuTime = 0;

    DMatrix dm = DMatrix::Random(4000,4000);
    SMatrix sm = denseToSparseBlock(dm, 0, 0, 4000, 4000);
    dm.resize(0,0);
    MySelfAdjointEigenSolver<DMatrix> es;
    
    startCPUTime = std::clock();
    es.computeBlock(sm, 0, 0, 4000, 4000);
    endCPUTime = std::clock();
    cpuTime = static_cast<double>(endCPUTime-startCPUTime)/CLOCKS_PER_SEC;
    std::cerr << "Time to diagonalize 4000x4000 matrix: " << cpuTime << "s" << std::endl;
#endif
}
