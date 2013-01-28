#ifndef DEC_H
#define DEC_H
#include "simplicialComplex.hpp"
#include <utility>
#include <assert.h>
#include <type_traits>
template <typename SC>
class DEC;

enum FormType {NO_FORM=0, PRIMAL_FORM=1, DUAL_FORM=2, BOTH_FORM=3};
template <int Dim_, FormType TypeIn_, int NIn_, FormType TypeOut_, int NOut_, bool isVector_=false>
struct form_traits{
    const static int Dim = Dim_;
    const static FormType TypeIn = TypeIn_;//Typein = -1 means that this should resolve to a vector
    const static FormType TypeOut = TypeOut_;
    const static int NIn = NIn_;
    const static int NOut = NOut_;
    const static bool isVector = isVector_;

};


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

    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
        Expression expr;

};


template <int Dim, typename DynamicVector,FormType Type1, int N1>
//Type1 doesn't do anything
class Form: public FormExpression<form_traits<Dim,Type1,-1,Type1,N1,true>, DynamicVector>
{
    private:
        typedef Form<Dim,DynamicVector,Type1,N1> MyType;
    public:
        typedef form_traits<Dim,Type1,-1,Type1,N1,true> Traits;
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
            typedef FormExpression<Traits,DynamicVector> Parent;
        using Parent::expr;
        //typedef typename Parent::ExpressionType ExpressionType;
        Form(int size=0): Parent(DynamicVector::Zero(size))/*, expr(size)*/ {}
        Form(const DynamicVector & other): Parent(expr)
    {
        assert(expr.size() == other.size());
        expr = other;
    }
        Form(const MyType & other): Parent(expr)/*, expr(other.constData())*/ {}
        template <typename Traits1, typename Expr2>
            Form(FormExpression<Traits1,Expr2> const & rhs): Parent(rhs.expr)//, expr(rhs.expr)
        {
            static_assert(
                    (Traits1::NIn== -1) &&
                    (Traits1::TypeOut == Traits::TypeOut) &&
                    (Traits1::NOut == Traits::NOut)
                    , "Equality can't combine different forms");
            std::cout << rhs.expr.rows() << " " << rhs.expr.cols() << std::endl;
            std::cout << expr.rows() << " " << expr.cols() << std::endl;
            //std::cout << "Data: " << expr.transpose() << std::endl;
        }
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
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW
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

template <typename SimplicialComplex>
class FormFactory{
    protected:
        FormFactory(const SimplicialComplex & sc): m_sc(sc) {}

        template <FormType Type, int N>
            Form<SimplicialComplex::Dim, typename SimplicialComplex::NumTraits::DynamicVector,Type,N> genForm()
            {
                static_assert(N <= SimplicialComplex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
                return Form<SimplicialComplex::Dim,typename SimplicialComplex::NumTraits::DynamicVector,Type,N>(
                        (Type == PRIMAL_FORM) ?  m_sc.template numSimplices<N>() : m_sc.template numSimplices <SimplicialComplex::Dim-N>()
                        );
            }

    private:
        const SimplicialComplex & m_sc;
};

struct primal_tag{};
struct dual_tag{};

//TmD means Top dimension - D, in reality we're dealing with the Dth dimension 
//This foolery is to make sure that we are building operators by going up in dimension
//which is nice from a theoretical perspective (as the deRham complex goes upward)
template <typename SC, int TmD>
class HiddenOperatorContainer: public HiddenOperatorContainer<SC,TmD-1>
{

    friend DEC<SC>;
    protected:
        const static int TopD = SC::Dim;
        const static int D = TopD-TmD;
        typedef typename SC::NumTraits NumTraits;
        typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
        typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
        typedef HiddenOperatorContainer<SC,TmD-1> Parent;
        HiddenOperatorContainer(const SC & sc)
            : Parent(sc)
              , m_d_primal(sc.template b<D+1>().transpose())
              , m_d_dual(sc.template b<TmD>())
              , m_hodge_primal(sc.template numSimplices<D>())
              , m_hodge_dual(sc.template numSimplices<D>())
    {
        for(auto&& s: sc.template constSimplices<D>())
        {
            m_hodge_primal.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
            m_hodge_dual.data().diagonal()(s.Index()) = ((D*(TmD)%2==0)?1:-1)*s.Volume() / s.DualVolume();
        }

    }

    protected://Data
        FormOperator<form_traits<TopD,PRIMAL_FORM,D,PRIMAL_FORM,D+1>,SparseMatrixColMajor> m_d_primal;
        FormOperator<form_traits<TopD,DUAL_FORM,TmD,DUAL_FORM,TmD+1>,const SparseMatrixColMajor &> m_d_dual;

        //SparseMatrixColMajor m_d;
        FormOperator<form_traits<TopD,PRIMAL_FORM,D,DUAL_FORM,TmD>, DiagonalMatrix> m_hodge_primal;//primal D -> dual TmD
        FormOperator<form_traits<TopD,DUAL_FORM,TmD,PRIMAL_FORM,D>, DiagonalMatrix> m_hodge_dual;//dual TmD -> primal D

        /*
           template <FormType Type>
           const typename std::enable_if<Type == PRIMAL_FORM, decltype(m_hodge_primal)>::type &
           internal_h()
           {
           return m_hodge_primal;
           }
           */

    private:
        auto internal_h(primal_tag) const
            -> const decltype(m_hodge_primal) &
            {
                return m_hodge_primal;
            }
        auto  internal_h(dual_tag) const
            -> const decltype(m_hodge_dual) &
            {
                return m_d_dual;
            }
        auto internal_d(primal_tag) const
            -> const decltype(m_d_primal) &
            {
                return m_d_primal;
            }
        auto  internal_d(dual_tag) const
            -> const decltype(m_d_dual) &
            {
                return m_d_dual;
            }
    protected:
        template <FormType Type>
            auto internal_d() const
            -> const typename std::conditional<Type==PRIMAL_FORM, decltype(m_d_primal), decltype(m_d_dual) >::type &
            {
                return internal_d(typename std::conditional<Type==PRIMAL_FORM, primal_tag, dual_tag>::type());
            }
        template <FormType Type>
            auto internal_h() const
            -> const typename std::conditional<Type==PRIMAL_FORM, decltype(m_hodge_primal), decltype(m_hodge_dual) >::type &
            {
                return internal_h(typename std::conditional<Type==PRIMAL_FORM, primal_tag, dual_tag>::type());
            }
};
/*
   template <typename SC, int TmD>
   auto HiddenOperatorContainer<SC,TmD>::h_()
   -> const decltype(m_hodge_dual) & 
   {
   return m_hodge_dual;
   }
   */

template <typename SC>
class HiddenOperatorContainer<SC,0>{
    friend DEC<SC>;
    typedef typename SC::NumTraits NumTraits;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    protected:
    HiddenOperatorContainer(const SC & sc)
        : m_d_dual(sc.template b<SC::Dim>())
          , m_hodge_primal(sc.template numSimplices<SC::Dim>())
          , m_hodge_dual(sc.template numSimplices<SC::Dim>())
    {
        for(auto&& s: sc.template constSimplices<SC::Dim>())
        {
            m_hodge_primal.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
            m_hodge_dual.data().diagonal()(s.Index()) = s.Volume() / s.DualVolume();
        }
    }
    protected:

    protected:
    bool m_d_primal=false;
    FormOperator<form_traits<SC::Dim,DUAL_FORM,0,DUAL_FORM,1>,const SparseMatrixColMajor &> m_d_dual;
    FormOperator<form_traits<SC::Dim,PRIMAL_FORM,0,DUAL_FORM,SC::Dim>, DiagonalMatrix> m_hodge_primal;//primal D -> dual TmD
    FormOperator<form_traits<SC::Dim,DUAL_FORM,SC::Dim,PRIMAL_FORM,0>, DiagonalMatrix> m_hodge_dual;//dual TmD -> primal D
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
        auto internal_d(primal_tag) const//this will always fail if you do anything with the returned operator... should change it to an assert later though
            -> const decltype(m_d_primal) &
            {
                return m_d_primal;
            }
        auto  internal_d(dual_tag) const
            -> const decltype(m_d_dual) &
            {
                return m_d_dual;
            }
    protected:
        template <FormType Type>
            auto internal_d() const
            -> const typename std::conditional<Type==PRIMAL_FORM, decltype(m_d_primal), decltype(m_d_dual) >::type &
            {
                return internal_d(typename std::conditional<Type==PRIMAL_FORM, primal_tag, dual_tag>::type());
            }
    template <FormType Type>
        auto internal_h() const
        -> const typename std::conditional<Type==PRIMAL_FORM, decltype(m_hodge_primal), decltype(m_hodge_dual) >::type &
        {
            return internal_h(typename std::conditional<Type==PRIMAL_FORM, primal_tag, dual_tag>::type());
        }
};

template <typename SC>
class OperatorContainer: public HiddenOperatorContainer<SC,SC::Dim> {
    protected:
        OperatorContainer(const SC & sc)
            : HiddenOperatorContainer<SC,SC::Dim>(sc)
        {
        }
        typename SC::NumTraits::SparseMatrix m_hdhd;

};

template <typename SC>
class DEC: public FormFactory<SC>, public OperatorContainer<SC>
{
    public:
        typedef SC SimplicialComplex;
        typedef typename SC::NumTraits NumTraits;
        typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
        typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;

        static const int Dim = SC::Dim;

        typedef FormFactory<SC> FF;
        typedef OperatorContainer<SC> OC;

        DEC(const SimplicialComplex & sc)
            : FF(sc)
              , OC(sc)
              , m_sc(sc)
    {}


    private:
    public:
        template <int N, FormType Form = PRIMAL_FORM>
            auto d()const
            -> decltype(dynamic_cast<const HiddenOperatorContainer<SC,(Form==PRIMAL_FORM)?Dim-N:N> *>(this)->template internal_d<Form>())
            {    return dynamic_cast<const HiddenOperatorContainer<SC,(Form==PRIMAL_FORM)?Dim-N:N> *>(this)->template internal_d<Form>();}



        template <int N, FormType Form = PRIMAL_FORM>
            auto h()const
            -> decltype(dynamic_cast<const HiddenOperatorContainer<SC,(Form==PRIMAL_FORM)?Dim-N:N> *>(this)->template internal_h<Form>())
            {    return dynamic_cast<const HiddenOperatorContainer<SC,(Form==PRIMAL_FORM)?Dim-N:N> *>(this)->template internal_h<Form>();}
        template <typename Traits, typename Expr>
            auto d(const FormExpression<Traits, Expr> & rhs)
            -> decltype (d<Traits::NOut, Traits::TypeOut>() * rhs)
            {
                return 
                    d<Traits::NOut, Traits::TypeOut>() * rhs;
            }
        template <typename Traits, typename Expr>
            auto h(const FormExpression<Traits, Expr> & rhs)
            -> decltype (h<Traits::NOut, Traits::TypeOut>() * rhs)
            {
                return 
                    h<Traits::NOut, Traits::TypeOut>() * rhs;
            }

        /*
           template <typename Traits, typename Expr>
           auto h(const FormExpression<Traits, Expr> & rhs)
           -> decltype (d<Traits::NOut, Traits::TypeOut>() * rhs)
           {
           return 
           h<Traits::NOut, Traits::TypeOut>() * rhs;
           }
           */

        /*
           template <FormType Type, int N>
           const Form<SC,Type,N+1> d(const Form<SC,Type,N> & f)const
           {
           return Form<SC,Type,N+1>(d<N>()*f);
           }
           template <FormType Type, int N, typename Expression>
           auto d(const FormObject & f) const
           -> FormExpression<Type,N+1,decltype(d<N>()*f)>
           {
           typedef FormExpression<
           DUAL_FORM, N+1,
           decltype(d<N>()*f)
           > ResultExprType;
        //should map PRIMAL,N to DUAL,Dim-N
        //or it maps DUAL, N to PRIMAL,Dim-N
        return ResultExprType(d<N>() * f);
        }
        */
        /*
           template <int N, typename Mat>
           auto h(const Mat & m) const -> decltype(h<N>() * m)
           {
           return h<N>() * m;
           }
           */
        /*
           template <int N, typename Expression>
           auto h(const FormExpression<PRIMAL_FORM,N,Expression> & f) const
           -> FormExpression<DUAL_FORM,Dim-N,decltype(h<N>()*f)>
           {
           typedef FormExpression<
           DUAL_FORM, Dim-N,
           decltype(h<N>()*f)
           > ResultExprType;
        //should map PRIMAL,N to DUAL,Dim-N
        //or it maps DUAL, N to PRIMAL,Dim-N
        return ResultExprType(h<N>() * f);
        }
        template <int N, typename Expression>
        auto h(const FormExpression<DUAL_FORM,N,Expression> & f) const
        -> FormExpression<PRIMAL_FORM,Dim-N,decltype(((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f)>
        {
        typedef FormExpression<
        PRIMAL_FORM,Dim-N,
        decltype(
        ((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f
        )
        > ResultExprType;

        //should map PRIMAL,N to DUAL,Dim-N
        //or it maps DUAL, N to PRIMAL,Dim-N
        return ResultExprType(((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f);
        }
        */
        /*
           template <int N>
           const Form<SC,DUAL_FORM,Dim-N> h(const Form<SC,PRIMAL_FORM,N> & f)const
           {
           return h<N>() * f;
           }

           template <int N>
           const Form<SC,PRIMAL_FORM,N> h(const Form<SC,DUAL_FORM,Dim-N> & f)const
           {
        //enforce inverse hodge keeeps  *^{-1}* = -1^{k*n-k}
        return ((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f;
        }
        */


    private:
        const SimplicialComplex & m_sc;
};
#include "dec_overloads.hpp"
#endif
