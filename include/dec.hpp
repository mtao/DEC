#ifndef DEC_H
#define DEC_H
#include "simplicialComplex.hpp"
#include <utility>
#include <assert.h>
#include <type_traits>
enum FormType {PRIMAL, DUAL};

template <FormType Type1, int N1, typename Expression>
struct FormExpression//: public Expression
{
    public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    typedef Expression ExpressionType;
    static const int N = N1;
    static const FormType Type = Type1;
    /*
    FormExpression(const Expression & other): Expression(other) {}
    FormExpression(const FormExpression & other): Expression(static_cast<const Expression>(other)) {}
    */
    private:
    const Expression & expr;

};


template <typename SimplicialComplex,FormType Type1, int N1>
class Form: public SimplicialComplex::NumTraits::DynamicVector
//class Form: public FormExpression<Type1,N1, typename SimplicialComplex::NumTraits::DynamicVector>
{
private:
    typedef Form<SimplicialComplex,Type1,N1> MyType;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static const int N = N1;
    static const FormType Type = Type1;
    typedef typename SimplicialComplex::NumTraits NumTraits;
    typedef typename NumTraits::Scalar Scalar;
    //typedef FormExpression<Type,N,typename SimplicialComplex::NumTraits::DynamicVector> Parent;
        typedef typename NumTraits::DynamicVector Parent;
    //typedef typename Parent::ExpressionType ExpressionType;
    Form(): Parent(0) {}
    Form(const SimplicialComplex & sc)
        : Parent(//static_cast<const Parent>(
                     //typename Parent::ExpressionType(
                         (Type == PRIMAL) ?  sc.template numSimplices<N>() : sc.template numSimplices <SimplicialComplex::Dim-N>()
                                             )
                     //)
                 //)
    {
        static_assert(N <= SimplicialComplex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
        this->setConstant(Scalar(0));
    }
    Form(const typename SimplicialComplex::NumTraits::DynamicVector & other)
    {
        assert(this->size() == other.size());
        dynamic_cast<MyType >(*this) = other;
    }
    Form(const MyType & other): Parent(static_cast<const Parent>(other)) {}
    void init(const SimplicialComplex & sc)
    {
        resize((Type == PRIMAL) ?  sc.template numSimplices<N>() :sc.template numSimplices<SimplicialComplex::Dim-N>()
                                   );
        setConstant(Scalar(0));

    }
    template <FormType Type2, int N2, typename Expr2>
    MyType & operator=(const FormExpression<Type2,N2,Expr2> & rhs)
    {
        static_assert(Type == Type2, "Equality can't combine different forms");
        static_assert(N == N2, "Equality can't combine different forms");
        std::cout << this << std::endl;
        std::cout << this->transpose() << std::endl;
        static_cast<typename Parent::ExpressionType>(*this) = rhs.expr;
        return *this;
    }
};

template <typename SimplicialComplex>
class FormFactory{
protected:
    FormFactory(const SimplicialComplex & sc): m_sc(sc) {}

    template <FormType Type, int N>
    Form<SimplicialComplex,Type,N> genForm()
    {
        return Form<SimplicialComplex,Type,N>(m_sc);
    }

private:
    const SimplicialComplex & m_sc;
};


//TmD means Top dimension - D, in reality we're dealing with the Dth dimension 
//This foolery is to make sure that we are building operators by going up in dimension
//which is nice from a theoretical perspective (as the deRham complex goes upward)
template <typename SC, int TmD>
class HiddenOperatorContainer: public HiddenOperatorContainer<SC,TmD-1>
{

protected:
    const static int TopD = SC::Dim;
    const static int D = TopD-TmD;
    typedef typename SC::NumTraits NumTraits;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef HiddenOperatorContainer<SC,TmD-1> Parent;
    HiddenOperatorContainer(const SC & sc)
        : Parent(sc)
        , m_d(sc.template b<D+1>().transpose())
        , m_hodge(sc.template numSimplices<D>())
    {
        for(auto&& s: sc.template constSimplices<D>())
        {
            m_hodge.diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }

    }

protected://Data
    SparseMatrixColMajor m_d;
    DiagonalMatrix m_hodge;
};

template <typename SC>
class HiddenOperatorContainer<SC,0>{
protected:
    HiddenOperatorContainer(const SC &) {}
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

    static const unsigned int Dim = SC::Dim;

    typedef FormFactory<SC> FF;
    typedef OperatorContainer<SC> OC;

    DEC(const SimplicialComplex & sc)
        : FF(sc)
        , OC(sc)
        , m_sc(sc)
    {}


    template <int N>
    const SparseMatrixColMajor & d()const
    {
        return HiddenOperatorContainer<SC,Dim-N>::m_d;
    }

    template <int N, FormType Form = PRIMAL>
    const DiagonalMatrix & h()const
    {
        if(Form == PRIMAL)
        {
            return HiddenOperatorContainer<SC,Dim-N>::m_hodge;
        }
        else
        {
            return ((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse();
        }
    }

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
                DUAL, N+1,
                decltype(d<N>()*f)
                > ResultExprType;
        //should map Primal,N to Dual,Dim-N
        //or it maps Dual, N to Primal,Dim-N
        return ResultExprType(d<N>() * f);
    }
    /*
    template <int N, typename Mat>
    auto h(const Mat & m) const -> decltype(h<N>() * m)
    {
        return h<N>() * m;
    }
    */
    /*
    template <int N, typename Expression>
    auto h(const FormExpression<PRIMAL,N,Expression> & f) const
    -> FormExpression<DUAL,Dim-N,decltype(h<N>()*f)>
    {
        typedef FormExpression<
                DUAL, Dim-N,
                decltype(h<N>()*f)
                > ResultExprType;
        //should map Primal,N to Dual,Dim-N
        //or it maps Dual, N to Primal,Dim-N
        return ResultExprType(h<N>() * f);
    }
    template <int N, typename Expression>
    auto h(const FormExpression<DUAL,N,Expression> & f) const
    -> FormExpression<PRIMAL,Dim-N,decltype(((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f)>
    {
        typedef FormExpression<
                PRIMAL,Dim-N,
                decltype(
                    ((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f
                    )
                > ResultExprType;

        //should map Primal,N to Dual,Dim-N
        //or it maps Dual, N to Primal,Dim-N
        return ResultExprType(((N*(Dim-N)%2==0)?1:-1)*h<N>().inverse() * f);
    }
    */
    /*
        template <int N>
            const Form<SC,DUAL,Dim-N> h(const Form<SC,PRIMAL,N> & f)const
            {
                return h<N>() * f;
            }

        template <int N>
            const Form<SC,PRIMAL,N> h(const Form<SC,DUAL,Dim-N> & f)const
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
