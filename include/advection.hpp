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
    typedef typename DEC::Nm1Form VelocityFormType;
    typedef typename Complex::DimTraits DimTraits;
    typedef typename Complex::Vector Vector;
    typedef typename Vector::Scalar Scalar;
    typedef typename DEC::SparseMatrixColMajor SparseMatrix;
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    Particle(const DEC & dec, const Vector & p = Vector::Zero(), const Simplex * simplex=0)//, const Vector & v = Vector::Zero())
        : m_dec(&dec)
        , m_pos(p)
        , m_simplex(simplex)
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
    Simplex * simplex() {return m_simplex;}
    const Simplex * simplex() const {return m_simplex;}
    /*
    void v() {
        m_dec->simplex()->whitney
    }
    */
private:
    const DEC * m_dec;
    Vector m_pos = Vector::Zero();
    //Vector m_vel = Vector::Zero();
    const Simplex * m_simplex = NULL;
public:

    void advectInPlace(const VelocityFormType & form, Scalar dt) {
        if(m_simplex == NULL) {
            std::cout << "Simplex is lost, projecting" << std::endl;
            project();
        }
        std::cout << "Starting in simplex: " << m_simplex->Index() << "    ===================="<< std::endl;
        Vector velocity;
        bool reverse = dt < 0;
        if(reverse) dt = -dt;
        Scalar eps = dt * 0.001;
        while(dt >0) {
            m_dec->getVelocityInPlace(m_pos,*m_simplex, form,velocity);
            if(reverse) velocity *= -1;
            auto&& basis = m_dec->complex().whitneyBasis(*m_simplex);
            int sind = 0;
            int traversed_edge = -1;
            Scalar timeToIntersection = std::numeric_limits<Scalar>::max();
            for(typename SparseMatrix::InnerIterator it(m_dec->complex().template b<DimTraits::Top>(), m_simplex->Index()); it; ++it, ++sind) {
                auto&& lower = m_dec->complex().template simplex<DimTraits::Top-1>(it.row());
                Scalar localtti = (lower.Center()-m_pos).dot(basis.col(sind)) / velocity.dot(basis.col(sind));
                if(localtti > eps  && localtti < timeToIntersection && localtti < dt ) {
                    //std::cout << "Update: " << localtti << " " << m_simplex->Index() << "-->" <<it.row() << std::endl;
                    traversed_edge = lower.Index();
                    timeToIntersection = localtti;
                }
            }
            std::cout << "Current time to intersection: " << timeToIntersection << std::endl;
            if(traversed_edge == -1) {
                m_pos += dt * velocity;
                mtao::projectToSimplexInPlace(m_dec->complex(),m_pos,*m_simplex);
                return;
            } else {
                std::cout <<"Chosen edge: " << traversed_edge << std::endl;
                std::cout << dt << "<--" << timeToIntersection<< std::endl;
                dt -= timeToIntersection;
                m_pos += timeToIntersection * velocity;
                typename SparseMatrix::InnerIterator it2(m_dec->template d<DimTraits::Top-1>().expr, traversed_edge);
                if(it2.row() == m_simplex->Index()) {
                    std::cout << "I see myself: " << it2.row() << std::endl;
                    ++it2;
                    std::cout << "I changed to: " << it2.row() << std::endl;
                }
                if(it2) {
                    m_simplex = &(m_dec->complex().simplex(it2.row()));
                    std::cout << *m_simplex << std::endl;
                    std::cout << "  Simplex changed to: " << m_simplex->Index() << std::endl;
                    mtao::projectToSimplexInPlace(m_dec->complex(),m_pos,*m_simplex);
                    if(dt < 0) {return;}
                    traversed_edge = true;
                } else {
                    std::cout << "WTF i just intersected with a boundary" << std::endl;
                    return;
                }
            }

        }//while loop

    }//function
    typename DEC::DECTraits::template form<PRIMAL_FORM,Simplex::Dim>::type activeSimplex() {
        auto form = m_dec->template genForm<PRIMAL_FORM, Simplex::Dim>();
        form.expr.setZero();
        form.expr(m_simplex->Index()) = 1;
        return form;
    }
};


/*
template <typename DEC>
class AdvectionDEC: public DEC
{
public:
    typedef typename DEC::Complex Complex;
    typedef typename Complex::DimTraits DimTraits;
    typedef typename DEC::Nm1Form VelocityFormType;
    typedef typename Complex::Vector Vector;
    typedef typename Vector::Scalar Scalar;

protected:




};
*/


/*
                    */
#endif
