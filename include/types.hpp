#ifndef TYPES_UTIL_H
#define TYPES_UTIL_H

//This is to make sure that forms dont change size by accident
//#define EIGEN_NO_AUTOMATIC_RESIZING
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include "index.hpp"
typedef unsigned int uint;


template <typename T, int DIM=Eigen::Dynamic>
struct NumericalTraits
{
    typedef T Scalar;
    static const int Dim = DIM;
    typedef mtao::IndexSet<Dim> IndexSet;
    typedef Eigen::Matrix<T,DIM,1> Vector;
    typedef Eigen::Matrix<T,3,1> Vector3;
    typedef Eigen::Matrix<T,2,1> Vector2;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> DynamicVector;

    typedef Eigen::Ref< Vector        > VectorRef;
    typedef Eigen::Ref< Vector3       > Vector3Ref;
    typedef Eigen::Ref< Vector2       > Vector2Ref;
    typedef Eigen::Ref< DynamicVector > DynamicVectorRef;

    typedef Eigen::Ref< const Vector        > VectorConstRef;
    typedef Eigen::Ref< const Vector3       > Vector3ConstRef;
    typedef Eigen::Ref< const Vector2       > Vector2ConstRef;
    typedef Eigen::Ref< const DynamicVector > DynamicVectorConstRef;

    //Assume that I'm dealing with square matrices or dynamic matrices
    typedef Eigen::Matrix<T,DIM,DIM> Matrix;
    typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> DynamicMatrix;
    typedef Eigen::DiagonalMatrix<T,Eigen::Dynamic> DiagonalMatrix;
    typedef Eigen::SparseMatrix<T> SparseMatrix;
    typedef Eigen::SparseMatrix<T, Eigen::RowMajor> SparseMatrixRowMajor;
    typedef Eigen::SparseMatrix<T, Eigen::ColMajor> SparseMatrixColMajor;
    typedef Eigen::Triplet<Scalar> Triplet;
};

typedef NumericalTraits<float,Eigen::Dynamic> NumericalTraitsXf;
typedef NumericalTraits<double,Eigen::Dynamic> NumericalTraitsXd;

#endif
