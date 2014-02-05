#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE eigenHelpers test
#include <boost/test/unit_test.hpp>

#include <ctime>
#include <iostream>
#include "eigenHelpers.hpp"

using namespace EigenHelpers;

SMatrix constructSparseMatrix ()
{
    SparseMatrix<double> sm(5,5);
    sm.reserve(VectorXi::Constant(5,5));
    int k=1;
    for (int i=0; i<5; i++)
        for (int j=0; j<5; j++)
        {
            sm.insert(i,j) = k;
            k++;
        }
    sm.makeCompressed();
    return sm;
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
