#ifndef DEC_H
#define DEC_H
#include "simplicialComplex.hpp"
#include <utility>
#include <assert.h>
#include <type_traits>
enum FormType {PRIMAL_FORM, DUAL_FORM, BOTH_FORM};

template <FormType TypeIn_, int NIn_, FormType TypeOut_, int NOut_, typename Expression>
struct FormExpression
{
    const static FormType TypeIn = TypeIn_;//Typein = -1 means that this should resolve to a vector
    const static FormType TypeOut = TypeOut_;
    const static int NIn = NIn_;
    const static int NOut = NOut_;
    FormExpression(const Expression & exp): expr(exp) {}

    const Expression & expr;

};


template <typename DynamicVector,FormType Type1, int N1>
//Type1 doesn't do anything
class Form: public FormExpression<Type1,-1,Type1,N1, DynamicVector>
{
private:
    typedef Form<DynamicVector,Type1,N1> MyType;
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
    static const int N = N1;
    static const FormType Type = Type1;
    typedef FormExpression<Type,-1,Type,N,DynamicVector> Parent;
    //typedef typename Parent::ExpressionType ExpressionType;
    Form(int size=0): Parent(m_data), m_data(size) {m_data.setConstant(typename DynamicVector::Scalar(0));}
    Form(const DynamicVector & other): Parent(m_data)
    {
        assert(m_data.size() == other.size());
        m_data = other;
    }
    Form(const MyType & other): Parent(m_data), m_data(other.constData()) {}
    template <typename Expr>
    Form(const Expr & other): Parent(m_data), m_data(other) {}
    template <FormType TypeIn, int NIn, FormType TypeOut, int NOut, typename Expr2>
        MyType & operator=(const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr2> & rhs)
        {
            static_assert(
                    (NIn== -1) &&
                    (TypeOut == Type) &&
                    (NOut == N)
                    , "Equality can't combine different forms");
            std::cout << m_data.transpose() << std::endl;
            m_data = rhs.expr;
            return *this;
        }
    DynamicVector & data(){return m_data;}
    const DynamicVector & constData() const {return m_data;}
private:
    DynamicVector m_data;
};



template <FormType TypeIn_, int NIn_, FormType TypeOut_, int NOut_, typename MatrixType>
class FormOperator: public FormExpression<TypeIn_,NIn_,TypeOut_,NOut_, MatrixType>
{
private:
    typedef FormOperator<TypeIn_, NIn_, TypeOut_,NOut_,MatrixType> MyType;
    public:
        EIGEN_MAKE_ALIGNED_OPERATOR_NEW;
        static const int NIn = NIn_;
        static const FormType TypeIn = TypeIn_;
        static const int NOut = NOut_;
        static const FormType TypeOut = TypeOut_;
        typedef  FormExpression<TypeIn_,NIn_,TypeOut_,NOut_, MatrixType> Parent;
    template <typename Expr>
    FormOperator(const Expr & other): Parent(m_data), m_data(other) {}
    FormOperator(const MyType & other): Parent(m_data), m_data(other.constData()) {}
        FormOperator(): Parent(m_data) {}
    template <FormType _TypeIn, int _NIn, FormType _TypeOut, int _NOut, typename Expr2>
        MyType & operator=(const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr2> & rhs)
        {
            static_assert(
                    (TypeIn == _TypeIn) &&
                    (NIn==_NIn) &&
                    (TypeOut == _TypeOut) &&
                    (NOut == _NOut)
                    , "Equality can't combine different forms");
            m_data = rhs.expr;
            return *this;
        }
    MatrixType & data(){return m_data;}
    const MatrixType & constData() const {return m_data;}
    private:
        MatrixType m_data;

};

template <typename SimplicialComplex>
class FormFactory{
    protected:
        FormFactory(const SimplicialComplex & sc): m_sc(sc) {}

        template <FormType Type, int N>
            Form<typename SimplicialComplex::NumTraits::DynamicVector,Type,N> genForm()
            {
        static_assert(N <= SimplicialComplex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
                return Form<typename SimplicialComplex::NumTraits::DynamicVector,Type,N>(
          (Type == PRIMAL_FORM) ?  m_sc.template numSimplices<N>() : m_sc.template numSimplices <SimplicialComplex::Dim-N>()
              );
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
            m_hodge.data().diagonal()(s.Index()) = s.DualVolume() / s.Volume();
        }

    }

protected://Data
    FormOperator<BOTH_FORM,D,BOTH_FORM,D+1,SparseMatrixColMajor> m_d;
    FormOperator<PRIMAL_FORM,D,DUAL_FORM,TmD, DiagonalMatrix> m_hodge;
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
    const FormOperator<BOTH_FORM,N,BOTH_FORM,N+1,SparseMatrixColMajor> & d()const
    {
        return HiddenOperatorContainer<SC,Dim-N>::m_d;
    }

    template <int N, FormType Form = PRIMAL_FORM>
    const FormOperator<BOTH_FORM,N,BOTH_FORM,Dim-N,DiagonalMatrix> & h()const
    {
        if(Form == PRIMAL_FORM)
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
                DUAL_FORM, N+1,
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
    auto h(const FormExpression<PRIMAL_FORM,N,Expression> & f) const
    -> FormExpression<DUAL_FORM,Dim-N,decltype(h<N>()*f)>
    {
        typedef FormExpression<
                DUAL_FORM, Dim-N,
                decltype(h<N>()*f)
                > ResultExprType;
        //should map Primal,N to Dual,Dim-N
        //or it maps Dual, N to Primal,Dim-N
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

        //should map Primal,N to Dual,Dim-N
        //or it maps Dual, N to Primal,Dim-N
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
