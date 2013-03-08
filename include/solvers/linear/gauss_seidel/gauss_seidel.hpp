#ifndef GAUSS_SEIDEL_H
#define GAUSS_SEIDEL_H
#include "../linear.h"

/*
 * Ax=b
 */
template <typename MatrixType, typename VectorType>
struct DenseGaussSeidelCapsule: public LinearSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef LinearSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::A ;
    using CapsuleBase::b ;
    using CapsuleBase::x ;
    DenseGaussSeidelCapsule(const Matrix & A, const Vector & b, Vector & x):
        CapsuleBase(A,b,x) {}
    void step()
    {
        x = (A.template triangularView<Eigen::Lower>().solve(
                    b-A.template triangularView<Eigen::StrictlyUpper>()*x));
    }

};

template <typename Matrix, typename Vector>
void DenseGaussSeidelSolve(const Matrix & A, const Vector & b, Vector & x)
{
    auto solver = IterativeLinearSolver<DenseGaussSeidelCapsule<Matrix,Vector> >(A.rows(), 0.001);
    solver.solve(A,b,x);
}

template <typename MatrixType, typename VectorType>
struct SparseGaussSeidelCapsule: public LinearSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef LinearSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::A ;
    using CapsuleBase::b ;
    using CapsuleBase::x ;
    SparseGaussSeidelCapsule(const Matrix & A, const Vector & b, Vector & x):
        CapsuleBase(A,b,x)
    {
        invdiag.resize(x.rows());
        invdiag.setZero();
        for(int i=0; i < A.outerSize(); ++i)
            if(std::abs(A.coeff(i,i))>0.00001)
            invdiag(i)=1/A.coeff(i,i);
    }
    //Scalar error(){return err;}
    void step()
    {
        Scalar dx;
        for(int i=0; i< A.outerSize(); ++i)
        {
           dx = -b(i) + A.innerVector(i).dot(x);
           x(i) -= dx*invdiag(i);
           //err = std::max(std::abs(dx/diag(i)),err);
        }
    }
    private:
    Vector invdiag;
    Scalar err;

};
template <typename Matrix, typename Vector>
void SparseGaussSeidelSolve(const Matrix & A, const Vector & b, Vector & x)
{
    auto solver = IterativeLinearSolver<SparseGaussSeidelCapsule<Matrix,Vector> >(A.rows(), 0.001);
    solver.solve(A,b,x);
}
#endif
