#ifndef QUADRATIC_H
#define QUADRATIC_H
#include "geometry/util/types.h"
#include <memory>

/*
0 = Mz+q
min_z  z'Mz + q'z + c
  */
template <typename Matrix, typename Vector>
struct QuadraticSolverCapsule
{
    QuadraticSolverCapsule(const Matrix & M_, const Vector & q_, Vector & z_): M(M_), q(q_), z(z_) {}
    typedef typename Vector::Scalar Scalar;
    virtual void step() = 0;
    virtual Scalar error()
    {
        w = M*z-q;
        //std::cout << w.cwiseProduct(z).template lpNorm<Eigen::Infinity>() << " " << w.array().min(Scalar(0)).sum() << " " << z.array().min(Scalar(0)).sum() << std::endl;
        return w.cwiseProduct(z).template lpNorm<Eigen::Infinity>()
                - w.array().min(Scalar(0)).sum()
                - z.array().min(Scalar(0)).sum()
                ;
    }
    virtual void setTolerance(const Scalar eps){epsilon = eps;}

protected:
    const Matrix & M;
    const Vector & q;
    Vector w;
    Vector & z;
    Scalar epsilon;
};


template <typename DerivedCapsule, typename Matrix, typename Vector>
DerivedCapsule createQuadraticCapsule(const Matrix & M, const Vector & q, Vector & z)
{
    return DerivedCapsule(M,q,z);
}



template <typename Capsule>
struct IterativeQuadraticSolver
{
    typedef typename Capsule::Matrix Matrix;
    typedef typename Capsule::Vector Vector;
    typedef typename Vector::Scalar Scalar;
    IterativeQuadraticSolver(uint max_its=1000, Scalar eps=0.001):
        epsilon(eps), max_iterations(max_its) {}
    void solve()
    {
        if(!capsule){return;}
        uint iterations = 0;
        while(++iterations < max_iterations && capsule->error() > epsilon) {
            capsule->step();
        }
        if(capsule->error() > epsilon)
            std::cout << iterations << "/" << max_iterations << " " << capsule->error() << "/" << epsilon << std::endl;
    }
    void init(const Matrix & M, const Vector & q, Vector & z)
    {
        if(!capsule)
            capsule.reset(new Capsule(M,q,z));
        capsule->setTolerance(epsilon);
    }
    void solve(const Matrix & M, const Vector & q, Vector & z)
    {
        init(M,q,z);
        solve();
    }
    bool success()
    {return (epsilon >= capsule->error());}
    Capsule * getCapsule(){return capsule.get();}

private:
    std::unique_ptr<Capsule> capsule;
    Scalar epsilon;
    uint max_iterations;
    bool inited = false;
};
#endif // QUADRATIC_H
