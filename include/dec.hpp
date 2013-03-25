#ifndef _DEC_H_
#define _DEC_H_
#include "simplicialComplex.hpp"
#include <utility>
#include <assert.h>
#include <type_traits>
#include <array>
#include "types.hpp"
template <typename Complex_, bool Interior>
class DEC;
namespace mtao_internal {



template <typename Traits_, typename Expression_>
struct FormExpression
{
    typedef Traits_ Traits;
    typedef Expression_ Expression;
    const static int Dim = Traits::Dim;
    const static FormType TypeIn = Traits::TypeIn;//Typein = -1 means that this should resolve to a vector
    const static FormType TypeOut = Traits::TypeOut;
    const static int NIn = Traits::NIn;
    const static int NOut = Traits::NOut;
    typedef Expression ExpressionType;
    FormExpression(const Expression & exp): expr(exp) {}
    Expression & data(){return expr;}
    const Expression & constData() const {return expr;}

    Expression expr;

};


template <int Dim, typename DynamicVector,FormType Type1, int N1>
//Type1 doesn't do anything
class Form: public FormExpression<form_operator_traits<Dim,Type1,-1,Type1,N1,true>, DynamicVector>
{
private:
    typedef Form<Dim,DynamicVector,Type1,N1> MyType;
public:
    typedef form_operator_traits<Dim,Type1,-1,Type1,N1,true> Traits;
    typedef typename DynamicVector::Scalar Scalar;
    typedef FormExpression<Traits,DynamicVector> Parent;
    using Parent::expr;
    //typedef typename Parent::ExpressionType ExpressionType;
    Form(int size=0): Parent(DynamicVector::Zero(size))/*, expr(size)*/ {}
    Form(const DynamicVector & other): Parent(expr)
    {
        assert(expr.size() == other.size());
        expr = other;
    }
    Form(const MyType & other): Parent(expr)/*, expr(other.constData())*/ {
        this->expr = other.expr;
    }
    template <typename Traits1, typename Expr2>
    Form(FormExpression<Traits1,Expr2> const & rhs): Parent(rhs.expr)//, expr(rhs.expr)
    {
        static_assert(
                    (Traits1::NIn== -1) &&
                    (Traits1::TypeOut == Traits::TypeOut) &&
                    (Traits1::NOut == Traits::NOut)
                    , "Equality can't combine different forms");
    }
    std::vector<std::array<float, N1+1> > toRenderableData() {
        std::vector<std::array<float, N1+1> > ret(expr.rows());
        for(int i=0; i < expr.rows(); ++i) {
            ret[i].fill(static_cast<float>(expr[i]));
        }
        return ret;
    }
    Scalar & operator()(int i){return expr(i);}
    const Scalar & operator()(int i) const {return expr(i);}

private:
    //DynamicVector expr;
};



template <typename Traits_, typename MatrixType>
class FormOperator: public FormExpression<Traits_, MatrixType>
{
private:
    typedef FormOperator<Traits_,MatrixType> MyType;
public:
    typedef Traits_ Traits;
    typedef  FormExpression<Traits, MatrixType> Parent;
    using Parent::expr;
    //template <typename Expr>
    //FormOperator(const Expr & other): Parent(expr){}//, expr(other) {}
    FormOperator(const MyType & other): Parent(other.expr)/*, expr(other.constData())*/ {}
    //FormOperator(): Parent(expr) {}
    FormOperator(const MatrixType & mat): Parent(mat) {}
    template <typename _Traits, typename Expr2>
    FormOperator(const FormExpression<_Traits,Expr2> & rhs): Parent(expr)
    {
        static_assert(
                    (Traits::TypeIn == _Traits::TypeIn) &&
                    (Traits::NIn==_Traits::NIn) &&
                    (Traits::TypeOut == _Traits::TypeOut) &&
                    (Traits::NOut == _Traits::NOut)
                    , "Equality can't combine different forms");
        expr = rhs.expr;
    }
private:
    //MatrixType m_data;

};

template <typename DECTraits>
class FormFactory{
protected:
    typedef typename DECTraits::NumTraits NumTraits;
    typedef typename DECTraits::DimTraits DimTraits;
    typedef typename DECTraits::Complex Complex;
    FormFactory(const Complex & sc): m_sc(sc) {}

public:
    template <FormType Type = FormType::Primal, int N = Complex::Dim>
    Form<DimTraits::Dim, typename Complex::NumTraits::DynamicVector,Type,N> genForm() const
    {
        static_assert(N <= Complex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
        return Form<Complex::Dim,typename Complex::NumTraits::DynamicVector,Type,N>(
                    (Type == FormType::Primal) ?  m_sc.template numSimplices<N>() : m_sc.template numSimplices <Complex::Dim-N>()
                                             );
    }

private:
    const Complex & m_sc;
};

struct primal_tag{};
struct dual_tag{};
struct dual_interior_tag{};

template <typename DECTraits>
class OperatorContainerPrivateBase {
    friend class DEC<typename DECTraits::Complex, DECTraits::Interior>;
    typedef typename DECTraits::template operator_private<DECTraits::TopDim> MyTraits;

    public:
    typedef typename DECTraits::NumTraits NumTraits;
    typedef typename DECTraits::DimTraits DimTraits;
    typedef typename DECTraits::Complex Complex;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
protected:
    OperatorContainerPrivateBase(const Complex & sc)
    // : m_d_dual(sc.template b<Complex_::Dim>())
        : m_hodge_primal(
              sc.template numSimplices<DimTraits::Dim>()
              )
        , m_hodge_dual(
              sc.template numSimplices<DimTraits::Top - DimTraits::Dim>()
              )
    {
        for(auto&& s: sc.template constSimplices<DimTraits::Dim>())
        {
            m_hodge_primal.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }
        for(auto&& s: sc.template constSimplices<DimTraits::Top-DimTraits::Dim>())
        {
            m_hodge_dual.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }
    }
protected:

protected:
    typename MyTraits::hodge_primal_type m_hodge_primal;
    typename MyTraits::hodge_dual_type m_hodge_dual;
private:
    auto internal_h(primal_tag) const
    -> const decltype(m_hodge_primal) &
    {
        return m_hodge_primal;
    }
    auto  internal_h(dual_tag) const
    -> const decltype(m_hodge_dual) &
    {
        return m_hodge_dual;
    }
public:
    template <FormType Type>
    auto internal_h() const
    -> const typename MyTraits::template hodge_type<Type>::type &
    {
        return internal_h(typename std::conditional<Type==FormType::Primal, primal_tag, dual_tag>::type());
    }
};




//TmD means Top dimension - D, in reality we're dealing with the Dth dimension
//This foolery is to make sure that we are building operators by going up in dimension
//which is nice from a theoretical perspective (as the deRham complex goes upward)
template <typename DECTraits, typename DimTraits>
class OperatorContainerPrivate: public DECTraits::template operator_private<DimTraits::Dim+1>::type
{

    friend class DEC<typename DECTraits::Complex, DECTraits::Interior>;
    typedef typename DECTraits::template operator_private<DimTraits::Dim> MyTraits;
protected:
    const static int TopD = DimTraits::Top;
    const static int D = DimTraits::Dim;
    const static int TmD = TopD - D;
    typedef typename DECTraits::NumTraits NumTraits;
    typedef typename DECTraits::Complex Complex;
    PLAIN_CLASS_NUM_DEFS
    typedef typename DECTraits::template operator_private<DimTraits::UpperTraits::Dim>::type Parent;
    OperatorContainerPrivate(const Complex & sc)
        : Parent(sc)
        , m_d_primal(sc.template b<D+1>().transpose())
        , m_d_dual(sc.template b<TmD>())
        , m_d_dual_interior(sc.template interior<TmD-1>() * sc.template b<TmD>())
        , m_hodge_primal(
              sc.template numSimplices<D>()
              )
        , m_hodge_dual(
              sc.template numSimplices<TmD>()
              )
    {
        for(auto&& s: sc.template constSimplices<D>())
        {
            m_hodge_primal.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }
        for(auto&& s: sc.template constSimplices<TmD>())
        {
            m_hodge_dual.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }

    }

protected://Data
    typename MyTraits::d_primal_type m_d_primal;
    typename MyTraits::d_dual_type m_d_dual;
    typename MyTraits::d_dual_interior_type m_d_dual_interior;

    typename MyTraits::hodge_primal_type m_hodge_primal;
    typename MyTraits::hodge_dual_type m_hodge_dual;


private:
    auto internal_h(primal_tag) const
    -> const decltype(m_hodge_primal) &
    {
        return m_hodge_primal;
    }
    auto  internal_h(dual_tag) const
    -> const decltype(m_hodge_dual) &
    {
        return m_hodge_dual;
    }
    auto internal_d(primal_tag) const
    -> const decltype(m_d_primal) &
    {
        return m_d_primal;
    }
    auto internal_d(dual_interior_tag) const
    -> const decltype(m_d_dual_interior) &
    {
        return m_d_dual_interior;
    }
    auto  internal_d(dual_tag) const
    -> const decltype(m_d_dual) &
    {
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wreturn-stack-address"
        return m_d_dual;
#pragma GCC diagnostic pop
    }
protected:
    template <FormType Type, bool Interior>
    auto internal_d() const
    -> const typename MyTraits::template d_type<Type,Interior>::type &
    {
        return internal_d(typename std::conditional<Type==FormType::Primal, primal_tag
                          , typename std::conditional<Interior, dual_interior_tag, dual_tag>::type >::type());
    }
    template <FormType Type>
    auto internal_h() const
    -> const typename MyTraits::template hodge_type<Type>::type &
    {
        return internal_h(typename std::conditional<Type==FormType::Primal, primal_tag, dual_tag>::type());
    }
};


template <typename DECTraits>
class OperatorContainer: public DECTraits::template operator_private<0>::type{
protected:
    OperatorContainer(const typename DECTraits::Complex & sc)
        : DECTraits::template operator_private<0>::type(sc)
    {
    }
    typename DECTraits::Complex::NumTraits::SparseMatrix m_hdhd;

};
};

template <typename Complex_, bool Interior_ = false>
class DEC: public mtao_internal::template dec_traits<Complex_,Interior_>::factory_type
        , public mtao_internal::template dec_traits<Complex_,Interior_>::operator_type
{
public:
    typedef Complex_ Complex;
    const static bool Interior = Interior_;
    typedef mtao_internal::template dec_traits<Complex,Interior> DECTraits;
    typedef typename DECTraits::NumTraits NumTraits;
    typedef typename DECTraits::DimTraits DimTraits;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;

    static const int Dim = DimTraits::Dim;


    typedef typename DECTraits::factory_type FF;
    typedef typename DECTraits::operator_type OC;

    DEC(const Complex & sc)
        : FF(sc)
        , OC(sc)
        , m_sc(sc)
    {
    }
    void init() {
        FF::init(m_sc);
        OC::init(m_sc);
    }


public:
    template <int M, FormType Form = FormType::Primal, bool Int=Interior>
    auto d()const
    -> const typename DECTraits::template dec_operator<M,Form>::d_type &
    { return DECTraits::template dec_operator<M,Form>::type::template internal_d<Form,Int>();}




    template <int M, FormType Form = FormType::Primal>
    auto h()const
    -> const typename DECTraits::template dec_operator<M,Form>::h_type &
    { return DECTraits::template dec_operator<M,Form>::type::template internal_h<Form>();}


    template <typename Traits, typename Expr>
    auto d(const mtao_internal::FormExpression<Traits, Expr> & rhs)
    -> decltype (d<Traits::NOut, Traits::TypeOut>() * rhs)
    {
        return
                d<Traits::NOut, Traits::TypeOut>() * rhs;
    }
    template <typename Traits, typename Expr>
    auto h(const mtao_internal::FormExpression<Traits, Expr> & rhs)
    -> decltype (h<Traits::NOut, Traits::TypeOut>() * rhs)
    {
        return
                h<Traits::NOut, Traits::TypeOut>() * rhs;
    }

    const Complex & complex() const {return m_sc;}
    //====================================================
    //====================================================
    //==               Form Generation                  ==
    //====================================================
    //====================================================
    template <FormType Type=FormType::Primal, int M=0>
    typename DECTraits::template form<Type,M>::type genForm() const {
        return FF::template genForm<Type,M>();
    }

    //====================================================
    //====================================================
    //==            Interpolation                       ==
    //====================================================
    //====================================================


    typedef typename DECTraits::template form<FormType::Primal, Dim-1>::type Nm1Form;
    typedef typename Complex::SCTraits::template internal_complex<Dim>::type::WhitneyBasis VelocityBasis;
    typedef typename Complex::SCTraits::template internal_complex<Dim>::type::WhitneyCoefficients VelocityCoefficients;
    typedef typename Complex::Vector Vector;
    void getVelocityInPlace(const Vector & p, const typename Complex::template TraitsContainer<Dim>::simplextype & simplex, const Nm1Form & form, Vector & v) const{
        VelocityCoefficients coeffs;
        auto&& b = m_sc.template b<Dim>();
        auto&& basis = m_sc.whitneyBasis(simplex);
        auto&& centers = m_sc.whitneyCenters(simplex);
        int i=0;
        for(typename SparseMatrixColMajor::InnerIterator it(b,simplex.Index()); it; ++it, ++i) {
            auto&& lower = m_sc.template simplex<Dim-1>(it.row());
            coeffs(i) =
                    it.value() *
                    form(lower.Index())
                    * (1-basis.col(i).dot(p-centers.col(i))/basis.col(i).squaredNorm());
        }
        v = m_sc.whitneyBasis(simplex) * coeffs;
    }
    Vector getVelocity(const Vector & p, const typename Complex::template TraitsContainer<Dim>::simplextype & simplex, const Nm1Form & form) const{
        Vector v;
        getVelocityInPlace(p,simplex,form,v);
        return v;

    }
    std::vector<Vector> velocityField(const Nm1Form & form) {
    std::vector<Vector> vfield(m_sc.numSimplices());
        for(auto&& s: m_sc.constSimplices()) {
            getVelocityInPlace(s.Center(), s,form,vfield[s.Index()]);
        }
        return vfield;
    }


private:
    const Complex & m_sc;
};
#include "dec_overloads.hpp"
#endif
