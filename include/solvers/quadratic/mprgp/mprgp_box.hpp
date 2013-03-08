#ifndef MPRGP_BOX_HPP
#define MPRGP_BOX_HPP

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
struct MPRGPBoxCapsule: public QuadraticSolverCapsule<MatrixType, VectorType>
{
    typedef MatrixType Matrix;
    typedef VectorType Vector;
    typedef typename Vector::Scalar Scalar;
    typedef QuadraticSolverCapsule<Matrix,Vector> CapsuleBase;
    using CapsuleBase::M ;
    using CapsuleBase::q ;
    using CapsuleBase::z ;
    using CapsuleBase::epsilon ;
    MPRGPBoxCapsule(const Matrix & M, const Vector & q, Vector & z):
        CapsuleBase(M,q,z)
    {
        int rows = q.rows();
        _theta = _beta = _theta_hat = _nu = Vector::Zero(rows);
        r = M*z-q;

        alpha_bar = AlphaBarFunction::f(M,epsilon);
    }
    void setL(const Vector & L)
    {
        l = L;
        active_set(z);
        p=_theta;
    }

    Scalar error()
    {
        return _nu.norm();
    }
    void step()
    {
        active_set(z);
        //std::cout << "vz: " <<z.transpose() << std::endl;
        //std::cout << "Errors: " << _beta.squaredNorm() << " " <<  _theta.dot(_theta_hat) << std::endl;
        if(_beta.squaredNorm() <= _theta.dot(_theta_hat))
        {
            Mp=M*p;
            alpha_cg = r.dot(p)/p.dot(Mp);
            alpha_f = std::numeric_limits<Scalar>::infinity();
            for(int i=0; i<z.rows(); ++i)
            {//max timestep in direction to not hit boundary
                if(p(i)>0)
                {//
                    alpha_f = std::min(-(z(i)-l(i))/p(i),alpha_f);
                }
                else if (p(i) < 0)
                {
                    alpha_f = std::min(-(l(i)-z(i))/p(i),alpha_f);
                }
            }
            //std::cout << alpha_cg << " " << alpha_f << std::endl;
            if(alpha_cg <= alpha_f)
            {//CG Step
                //std::cout << "CG\n";
                z -= alpha_cg * p;
                r -= alpha_cg * Mp;
                gamma = theta(z).dot(Mp)/p.dot(Mp);
                p=_theta-gamma*p;
                active_set(z);

            }
            else
            {//Expansion
                //std::cout << "EX\n";
                Scalar eval=z.dot(.5*M*z-q);
                z-=alpha_f*p;
                r-=alpha_f*Mp;
                z-=alpha_bar*theta(z);

                bindToL(z);
                if(eval < z.dot(.5*M*z-q))
                    alpha_bar/=2;
                r = M*z-q;
                active_set(z);
                p=_theta;

            }
        }
        else
        {

                //std::cout << "PR\n";
            Md = M*_beta;
            //std::cout << "theta: " << _theta.transpose() << std::endl;
            //std::cout << "Beta: " << _beta.transpose() << std::endl;
            alpha_cg = r.dot(_beta)/_beta.dot(Md);
            z-=alpha_cg * _beta;
            r -= alpha_cg * Md;
            active_set(z);
            p=_theta;

        }

        //std::cout << "^z: " <<z.transpose() << std::endl;

    }
private:
    Vector l;
    Scalar rdz, alpha_cg, alpha_f, gamma, alpha_bar;
    Vector p,Mp,r,_theta,_beta,_theta_hat, _nu,h,d,Md;
    inline void active_set(const Vector & x)
    {
        for(int i=0; i<x.rows(); ++i)
        {
            if(l(i) == 0)
            {
                _theta_hat(i) = 0;
                _theta(i) = 0;
                _beta(i) = 0;
                continue;
            }
            if(x(i)<l(i))
            {
                if(x(i) > -l(i))
                {//x is on interior
                    _theta(i)=r(i);
                    _beta(i) = 0;
                    if(r(i) > 0)
                    {
                        _theta_hat(i)=fmin((l(i)-x(i))/alpha_bar, _theta(i));
                    }
                    else
                    {
                        _theta_hat(i)=fmax(-(l(i)-x(i))/alpha_bar, _theta(i));
                    }
                }
                else
                {//x is on negative boudnary
                    _beta(i) = fmin(0,r(i));
                    _theta(i)=0;
                    _theta_hat(i)=0;
                }
            }
            else
            {//x is on positive boundary
                _beta(i) = fmax(0,r(i));
                _theta(i)=0;
                _theta_hat(i)=0;

            }
        }
        _nu = _theta + _beta;
    }
    inline void bindToL(Vector & x) {
        x.noalias() = (x.array() < 0).select(x.cwiseMax(-l),x.cwiseMin(l));
    }

    inline Vector & theta(const Vector & x)
    {
        for(int i=0; i<x.rows(); ++i)
            if(std::abs(x(i))<l(i))
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
            _beta(i) = (std::abs(x(i))>l(i))?0:r(i);
        _beta.noalias() = _beta.cwiseMin(0);
        return _beta;
    }
    inline Vector & nu(const Vector & x)
    {
        _nu = theta(x) + beta(x);
        return _nu;
    }

};

template <typename Matrix, typename Vector>
void SparseMPRGPBoxSolve(const Matrix & M, const Vector & q, Vector & z, const Vector & L)
{
    auto solver = IterativeQuadraticSolver<MPRGPBoxCapsule<Matrix,Vector,
            sparseAlphaBar<Matrix,typename Matrix::Scalar>
            > >(3*M.rows(), 0.001);

    /*
    solver.init(M,q,z);
    solver.getCapsule()->setL(L);
    solver.solve();
    if(!solver.success())
        */
    {
        z.setZero();
        solver.init(M,q,z);
        solver.getCapsule()->setL(L);
        solver.solve();
    }
}



#endif // MPRGP_BOX_HPP
