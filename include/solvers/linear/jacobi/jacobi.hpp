#ifndef JACOBI_H
#define JACOBI_H
#include <Eigen/Dense>


#include "../linear.h"

/*
 * Ax=b
 */


template <typename MatrixType, typename VectorType>
struct DenseJacobiCapsule: public LinearSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef LinearSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::A ;
    using CapsuleBase::b ;
    using CapsuleBase::x ;
    DenseJacobiCapsule(const Matrix & A, const Vector & b, Vector & x):
        CapsuleBase(A,b,x)
    {
        invdiag.resize(x.rows());
        for(int i=0; i<A.rows(); ++i)
        {
            invdiag(i) = A(i,i);
            if(std::abs(invdiag(i))>0.000001)
                invdiag(i)=1/invdiag(i);
        }
    }
    void step()
    {
        dx = (b-A*x).cwiseProduct(invdiag);
        x+=dx;
    }
    private:
    Vector invdiag, dx;

};

template <typename Matrix, typename Vector>
void DenseJacobiSolve(const Matrix & A, const Vector & b, Vector & x)
{
    auto solver = IterativeLinearSolver<DenseJacobiCapsule<Matrix,Vector> >(A.rows(), 0.001);
    solver.solve(A,b,x);
}




template <typename MatrixType, typename VectorType>
struct SparseJacobiCapsule: public LinearSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef LinearSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::A ;
    using CapsuleBase::b ;
    using CapsuleBase::x ;
    SparseJacobiCapsule(const Matrix & A, const Vector & b, Vector & x):
        CapsuleBase(A,b,x)
    {
        invdiag.resize(x.rows());
        for(int i=0; i<A.rows(); ++i)
        {
            invdiag(i) = A.coeff(i,i);
            if(std::abs(invdiag(i))>0.000001)
                invdiag(i)=1/invdiag(i);
        }
    }

    void step()
    {
        dx = (b-A*x).cwiseProduct(invdiag);
        x+=dx;
    }
    private:
    Vector invdiag, dx;
};

template <typename Matrix, typename Vector>
void SparseJacobiSolve(const Matrix & A, const Vector & b, Vector & x)
{
    auto solver = IterativeLinearSolver<SparseJacobiCapsule<Matrix,Vector> >(A.rows(), 0.001);
    solver.solve(A,b,x);
}


#endif
