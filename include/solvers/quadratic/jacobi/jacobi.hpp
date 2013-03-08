#ifndef JACOBI_HPP
#define JACOBI_HPP
#include "../quadratic.h"


/*
 * Ax=b
 */
template <typename MatrixType, typename VectorType>
struct DenseJacobiCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    DenseJacobiCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z) {}
    void step()
    {

        Scalar dz;
        z_ = M * z;
        for(int i=0; i<M.rows(); ++i)
        {
            if(M(i,i))
            {
                dz=-q(i)+z(i);
                z(i)-=dz/M(i,i);
                if(z(i)<0) z(i)=0;
            }
        }


    }
private:
    Vector z_;


};

template <typename Matrix, typename Vector>
void DenseJacobiSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<DenseJacobiCapsule<Matrix,Vector> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}

template <typename MatrixType, typename VectorType>
struct SparseJacobiCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    SparseJacobiCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z)
    {
        invdiag.resize(z.rows());
        invdiag.setZero();
        for(int i=0; i < M.outerSize(); ++i)
            if(std::abs(M.coeff(i,i))>0.00001)
                invdiag(i)=1/M.coeff(i,i);
    }
    void step()
    {
        z_ = M*z;
        Scalar dz;
        for(int i=0; i< M.outerSize(); ++i)
        {
            dz = -q(i) + z(i);
            z(i) -= dz*invdiag(i);
            if(z(i)<0) z(i) = 0;
        }
    }
private:
    Vector invdiag, z_;
    Scalar err;

};
template <typename Matrix, typename Vector>
void SparseJacobiSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<SparseJacobiCapsule<Matrix,Vector> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}


















template <typename MatrixType, typename VectorType, int MinimizedIterations, int FullIterations>
struct DenseJacobiSubspaceMinimizationCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    DenseJacobiSubspaceMinimizationCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z) {}
    void step()
    {

        Scalar dz;
        for(int k=0; k<FullIterations; ++k)
        {
            z_ = M * z;
            for(int i=0; i<M.rows(); ++i)
            {
                if(M(i,i))
                {
                    dz=-q(i)+z(i);
                    z(i)-=dz/M(i,i);
                    if(z(i)<0) z(i)=0;
                }
            }
        }
        for(int k=0; k<MinimizedIterations; ++k)
        {
            z_ = M * z;
            for(int i=0; i<M.rows(); ++i)
            {
                if(M(i,i))
                {
                    if(z(i) == 0) continue;
                    dz=-q(i)+z(i);
                    z(i)-=dz/M(i,i);
                    if(z(i)<0) z(i)=0;
                }
            }


        }
    }
private:
    Vector z_;

};

    template <typename Matrix, typename Vector>
void DenseJacobiSubspaceMinimizationSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<DenseJacobiSubspaceMinimizationCapsule<Matrix,Vector,3,5> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}

template <typename MatrixType, typename VectorType, int MinimizedIterations, int FullIterations>
struct SparseJacobiSubspaceMinimizationCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    SparseJacobiSubspaceMinimizationCapsule(const Matrix & M, const Vector & q, Vector & z):
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
        {
            z_ = M * z;
            for(int i=0; i< M.outerSize(); ++i)
            {
                dz = -q(i) + z(i);
                z(i) -= dz*invdiag(i);
                if(z(i)<0) z(i) = 0;
            }
        }
        for(int k=0; k < MinimizedIterations; ++k)
        {
            z_ = M * z;
            for(int i=0; i< M.outerSize(); ++i)
            {
                if(z(k)==0) continue;
                dz = -q(i) + z(i);
                z(i) -= dz*invdiag(i);
                if(z(i)<0) z(i) = 0;
            }
        }
     }
    private:
    Vector invdiag, z_;
    Scalar err;

};
    template <typename Matrix, typename Vector>
void SparseJacobiSubspaceMinimizationSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<SparseJacobiSubspaceMinimizationCapsule<Matrix,Vector,3,5> >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}











#endif // JACOBI_HPP
