#ifndef GAUSS_SEIDEL_HPP
#define GAUSS_SEIDEL_HPP
#include "../quadratic.h"


/*
 * Ax=b
 */
template <typename MatrixType, typename VectorType>
struct DenseGaussSeidelCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    DenseGaussSeidelCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z) {}
    void step()
    {

        Scalar dz;
        for(int i=0; i<M.rows(); ++i)
        {
            if(M(i,i))
            {
                dz=-q(i)+M.row(i).dot(z);
                z(i)-=dz/M(i,i);
                if(z(i)<0) z(i)=0;
            }
        }


    }

};

template <typename Matrix, typename Vector>
void DenseGaussSeidelSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<DenseGaussSeidelCapsule<Matrix,Vector> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}

template <typename MatrixType, typename VectorType>
struct SparseGaussSeidelCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    SparseGaussSeidelCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z)
    {
        invdiag.resize(z.rows());
        invdiag.setZero();
        for(int i=0; i < M.outerSize(); ++i)
            if(std::abs(M.coeff(i,i))>0.00001)
                invdiag(i)=1/M.coeff(i,i);
    }
    //Scalar error(){return err;}
    void step()
    {
        Scalar dz;
        for(int i=0; i< M.outerSize(); ++i)
        {
            dz = -q(i) + M.innerVector(i).dot(z);
            z(i) -= dz*invdiag(i);
            //err = std::max(std::abs(dx/diag(i)),err);
            if(z(i)<0) z(i) = 0;
        }
    }
private:
    Vector invdiag;
    Scalar err;

};
template <typename Matrix, typename Vector>
void SparseGaussSeidelSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<SparseGaussSeidelCapsule<Matrix,Vector> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}


















template <typename MatrixType, typename VectorType, int MinimizedIterations, int FullIterations>
struct DenseGaussSeidelSubspaceMinimizationCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    DenseGaussSeidelSubspaceMinimizationCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z) {}
    void step()
    {

        Scalar dz;
        for(int k=0; k<FullIterations; ++k)
            for(int i=0; i<M.rows(); ++i)
            {
                if(M(i,i))
                {
                    dz=-q(i)+M.row(i).dot(z);
                    z(i)-=dz/M(i,i);
                    if(z(i)<0) z(i)=0;
                }
            }
        for(int k=0; k<MinimizedIterations; ++k)
            for(int i=0; i<M.rows(); ++i)
            {
                if(M(i,i))
                {
                    if(z(i) == 0) continue;
                    dz=-q(i)+M.row(i).dot(z);
                    z(i)-=dz/M(i,i);
                    if(z(i)<0) z(i)=0;
                }
            }


    }

};

template <typename Matrix, typename Vector>
void DenseGaussSeidelSubspaceMinimizationSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<DenseGaussSeidelSubspaceMinimizationCapsule<Matrix,Vector,3,5> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}

template <typename MatrixType, typename VectorType, int MinimizedIterations, int FullIterations>
struct SparseGaussSeidelSubspaceMinimizationCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    SparseGaussSeidelSubspaceMinimizationCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z)
    {
        invdiag.resize(z.rows());
        invdiag.setZero();
        for(int i=0; i < M.outerSize(); ++i)
            if(std::abs(M.coeff(i,i))>0.00001)
                invdiag(i)=1/M.coeff(i,i);
    }
    //Scalar error(){return err;}
    void step()
    {
        Scalar dz;
        for(int k=0; k < FullIterations; ++k)
        for(int i=0; i< M.outerSize(); ++i)
        {
            dz = -q(i) + M.innerVector(i).dot(z);
            z(i) -= dz*invdiag(i);
            if(z(i)<0) z(i) = 0;
        }
        for(int k=0; k < MinimizedIterations; ++k)
        for(int i=0; i< M.outerSize(); ++i)
        {
                if(z(k)==0) continue;
            dz = -q(i) + M.innerVector(i).dot(z);
            z(i) -= dz*invdiag(i);
            if(z(i)<0) z(i) = 0;
        }
    }
private:
    Vector invdiag;
    Scalar err;

};
template <typename Matrix, typename Vector>
void SparseGaussSeidelSubspaceMinimizationSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<SparseGaussSeidelSubspaceMinimizationCapsule<Matrix,Vector,3,5> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}











#endif // GAUSS_SEIDEL_HPP
