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
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Particle(const DEC & dec, const Vector & p = Vector::Zero())//, const Vector & v = Vector::Zero())
        : m_dec(dec)
        , m_pos(p)
        //, m_vel(v)
    {
        //TODO: findNearestSimplex(vel);
        Scalar error = std::numeric_limits<Scalar>::max();
        auto&& sc = m_dec.complex();
        Vector offset = Vector::Zero();
        Eigen::Vector2f barycentric = Eigen::Vector2f(0,0);
//        std::cout << "Bary coords were: " << barycentric.transpose() << std::endl;
        for(auto&& simplex: sc.constSimplices()) {
            auto m = m_dec.complex().vertices(simplex);
            Eigen::Matrix<Scalar, Complex::EmbeddedDim, Complex::Dim-1> basis = m.rightCols(Complex::Dim-1);
            Vector origin = m.col(0);
            basis = basis - origin.rowwise().replicate(Complex::Dim-1);

            Eigen::Vector2f coords = m_dec.complex().barycentricCoords(simplex, m_pos);
            std::cout << coords.transpose() << std::endl;
            //REINSERT_HERE
            bool in_prism = true;
            for(int i=0; i < coords.rows(); ++i) {
                in_prism &= (coords(i) >= 0) && (coords(i) <= 1);
            }
            in_prism &= coords.sum() <= 1;

            if (in_prism) {// && Complex::EmbeddedDim == Complex::Dim+1) {
                Vector np = mtao::projectToSimplex(sc, p, simplex);
                Scalar err = (p-np).squaredNorm();
                if(err < error) {
                    error = err;
                    m_simplex = &simplex;
                    offset = np;
                    barycentric = coords;

                }
            }

        }
        std::cout << "Bary coords were: " << barycentric.transpose() << std::endl;
        std::cout << "Chosen offset: " << offset.transpose() << std::endl;
        m_pos.noalias() = offset;
        std::cout << "New position: " << m_pos.transpose() << std::endl;
    }
    Vector & p() {return m_pos;}
    const Vector & p() const {return m_pos;}
private:
    const DEC & m_dec;
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
