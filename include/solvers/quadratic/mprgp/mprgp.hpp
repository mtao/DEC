#ifndef MPRGP_HPP
#define MPRGP_HPP
#include "../quadratic.h"

template <typename MatrixType, typename Scalar>
struct denseAlphaBar
{
    static Scalar f(const MatrixType & M, Scalar epsilon)
    {
        Scalar ret = 0.0001;
        Scalar a;
        for(int j=0; j < M.rows(); ++j)
            for(int i=0; i < M.cols(); ++i)
            {
                a = std::abs(M(i,j));
                if(a > epsilon && a > ret)
                {
                    ret = a;
                }
            }
        return 2/ret;
    }
};
template <typename MatrixType, typename Scalar>
struct sparseAlphaBar
{
    static Scalar f(const MatrixType & M, Scalar epsilon)
    {
        Scalar ret = 0.0001;
        Scalar a;
        for(int k=0; k < M.outerSize(); ++k)
            for(typename MatrixType::InnerIterator it(M,k); it; ++it)
            {
                a = std::abs(it.value());
                if(a > epsilon && a > ret)
                {
                    ret = a;
                }
            }
        return 2/ret;
    }
};


/*
 * min x^T(Ax-b)
 * x \geq l
 */
template <typename MatrixType, typename VectorType, typename AlphaBarFunction>
struct MPRGPCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    using CapsuleBase::epsilon ;
    MPRGPCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z)
    {
        int rows = q.rows();
        l = _theta = _beta = _theta_hat = _gp = Vector::Zero(rows);
        r = M*z-q;
        active_set(z);
        p=_theta;
        alpha_bar = AlphaBarFunction::f(M,epsilon);
    }
    void setL(const Vector & L)
    {
        l = L;
    }

    Scalar error()
    {
        return _gp.norm();
    }
    void step()
    {
        active_set(z);
        if(_beta.squaredNorm() <= _theta.dot(_theta_hat))
        {
            Mp=M*p;
            alpha_cg = r.dot(p)/p.dot(Mp);
            alpha_f = std::numeric_limits<Scalar>::infinity();
            for(int i=0; i<z.rows(); ++i)
            {//max timestep in direction to not hit boundary
                if(p(i)>0)
                {
                    alpha_f = std::min((z(i)-l(i))/p(i),alpha_f);//find distance to boundary
                }
            }
            if(alpha_cg <= alpha_f)
            {//CG Step
                z -= alpha_cg * p;
                r -= alpha_cg * Mp;
                gamma = theta(z).dot(Mp)/p.dot(Mp);
                p=_theta-gamma*p;
                active_set(z);

            }
            else
            {//Expansion
                Scalar eval=z.dot(.5*M*z-q);
                z-=alpha_f*p;
                r-=alpha_f*Mp;
                z-=alpha_bar*theta(z);
                z=z.cwiseMax(0);//Bind z to boundary
                if(eval < z.dot(.5*M*z-q))
                    alpha_bar/=2;
                r = M*z-q;
                active_set(z);
                p=_theta;

            }
        }
        else
        {

            Md = M*_beta;
            alpha_cg = r.dot(_beta)/_beta.dot(Md);
            z-=alpha_cg * _beta;
            r -= alpha_cg * Md;
            active_set(z);
            p=_theta;

        }


    }
private:
    Vector l;
    Scalar rdz, alpha_cg, alpha_f, gamma, alpha_bar;
    Vector p,Mp,r,_theta,_beta,_theta_hat, _gp,h,d,Md;
    inline void active_set(const Vector & x)
    {
        for(int i=0; i<x.rows(); ++i)
            if(x(i)>l(i))//not on boundary
            {
                _theta(i)=r(i);
                _beta(i) = 0;
                _theta_hat(i)=fmin((x(i)-l(i))/alpha_bar, _theta(i));//distance to boundary

            }
            else//on the boundary
            {
                _beta(i) = fmin(0,r(i));
                _theta(i)=0;
                _theta_hat(i)=0;

            }
        _beta.noalias() = _beta.cwiseMin(0);//0 if it goes through boundary
        _gp = _theta + _beta;
    }

    inline Vector & theta(const Vector & x)
    {
        for(int i=0; i<x.rows(); ++i)
            if(x(i)>l(i))
            {
                _theta(i)=r(i);
            }
            else
            {
                _theta(i)=0;

            }
        return _theta;
    }

    inline Vector & beta(const Vector & x)
    {

        for(int i=0; i<x.rows(); ++i)
            _beta(i) = (x(i)<l(i))?0:r(i);
        _beta.noalias() = _beta.cwiseMin(0);
        return _beta;
    }
    inline Vector & gp(const Vector & x)
    {
        _gp = theta(x) + beta(x);
        return _gp;
    }

};

template <typename Matrix, typename Vector>
void DenseMPRGPSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<MPRGPCapsule<Matrix,Vector,
            denseAlphaBar<Matrix,typename Matrix::Scalar>
            > >(2*M.rows(), 0.001);
    solver.solve(M,q,z);
}
template <typename Matrix, typename Vector>
void SparseMPRGPSolve(const Matrix & M, const Vector & q, Vector & z)
{
    auto solver = IterativeQuadraticSolver<MPRGPCapsule<Matrix,Vector,
            sparseAlphaBar<Matrix,typename Matrix::Scalar>
            > >(3*M.rows(), 0.001);

    solver.init(M,q,z);
    solver.solve();
    if(!solver.success())
    {
        z.setZero();
        solver.solve(M,q,z);
    }
}



#endif // MPRGP_HPP
