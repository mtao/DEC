#ifndef ADVECTION_H
#define ADVECTION_H
#include <limits>
#include "util.hpp"
#include <iostream>
template <typename DEC>
class Particle {
public:
    typedef typename DEC::Complex Complex;
    typedef typename Complex::NSimplex Simplex;
    typedef typename Complex::DimTraits DimTraits;
    typedef typename Complex::Vector Vector;
    typedef typename Vector::Scalar Scalar;
    typedef typename DEC::SparseMatrixColMajor SparseMatrix;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Particle(const DEC & dec, const Vector & p = Vector::Zero())//, const Vector & v = Vector::Zero())
        : m_dec(&dec)
        , m_pos(p)
        //, m_vel(v)
    {
    }
        void project() {
        //TODO: findNearestSimplex(vel);
        Scalar error = std::numeric_limits<Scalar>::max();
        auto&& sc = m_dec->complex();
        Vector offset = Vector::Zero();

        for(auto&& simplex: sc.constSimplices()) {
            bool inside = true;
            for(typename SparseMatrix::InnerIterator it(m_dec->complex().b(), simplex.Index()); it; ++it) {
                auto&& s1 = m_dec->complex().template simplex<1>(it.row());
                Vector normal = (simplex.Center() - s1.Center()).normalized();
                if((m_pos- s1.Center()).dot(normal)<0) {
                    inside = false;
                    break;
                }
            }

            if (inside) {
                Vector np = mtao::projectToSimplex(sc, m_pos, simplex);
                Scalar err = (m_pos-simplex.Center()).squaredNorm();
                if(err < error) {
                    error = err;
                    m_simplex = &simplex;
                    offset = np;

                }
            }

        }
        std::cout << error << std::endl;
        m_pos.noalias() = offset;
        if(error > 0.01) {
            m_pos = Vector::Random();
            project();
        }
    }
    Vector & p() {return m_pos;}
    const Vector & p() const {return m_pos;}
    /*
    void v() {
        m_dec->simplex()->whitney
    }
    */
private:
    const DEC * m_dec;
    Vector m_pos = Vector::Zero();
    //Vector m_vel = Vector::Zero();
    const Simplex * m_simplex = 0;

};


template <typename DEC>
class AdvectionDEC: public DEC
{
public:
    typedef typename DEC::Complex Complex;
    typedef typename Complex::DimTraits DimTraits;
    typedef typename DEC::Nm1Form VelocityFormType;
    typedef typename Complex::Vector Vector;
    typedef typename Vector::Scalar Scalar;

    void advectInPlace(const Particle<DEC> & p, Scalar dt) {
        while(dt > 0) {


        }

    }
protected:




};


/*
                    */
#endif
