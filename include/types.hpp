#ifndef _TYPES_UTIL_H_
#define _TYPES_UTIL_H_

//This is to make sure that forms dont change size by accident
//#define EIGEN_NO_AUTOMATIC_RESIZING
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/LU>
#include "index.hpp"
typedef unsigned int uint;
#define MTAO_DEBUG_VECTOR(x) std::cerr << #x << " = " << x.transpose() << std::endl;

#define PLAIN_CLASS_NUM_DEFS \
    typedef typename NumTraits::Vector Vector;\
    typedef typename NumTraits::Scalar Scalar;\
    typedef typename NumTraits::Triplet Triplet;\
    typedef typename NumTraits::SparseMatrix SparseMatrix;\
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;\
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
namespace mtao_internal {
template <typename T, int DIM=Eigen::Dynamic>
struct NumTraits 
{
    typedef T Scalar;
    enum {dim = DIM};
    typedef mtao::IndexSet<Dim> IndexSet;
    typedef Eigen::Matrix<T,DIM,1> Vector;
    typedef Eigen::Matrix<T,3,1> Vector3;
    typedef Eigen::Matrix<T,2,1> Vector2;
    typedef Eigen::Matrix<T,Eigen::Dynamic,1> DynamicVector;

    typedef Eigen::Ref< Vector        > VectorRef;
    typedef Eigen::Ref< Vector3       > Vector3Ref;
    typedef Eigen::Ref< Vector2       > Vector2Ref;
    typedef Eigen::Ref< DynamicVector > DynamicVectorRef;

    typedef Eigen::Ref< const Vector        > VectorConstRef;
    typedef Eigen::Ref< const Vector3       > Vector3ConstRef;
    typedef Eigen::Ref< const Vector2       > Vector2ConstRef;
    typedef Eigen::Ref< const DynamicVector > DynamicVectorConstRef;

    //Assume that I'm dealing with square matrices or dynamic matrices
    typedef Eigen::Matrix<T,DIM,DIM> Matrix;
    typedef Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> DynamicMatrix;
    typedef Eigen::DiagonalMatrix<T,Eigen::Dynamic> DiagonalMatrix;
    typedef Eigen::SparseMatrix<T> SparseMatrix;
    typedef Eigen::SparseMatrix<T, Eigen::RowMajor> SparseMatrixRowMajor;
    typedef Eigen::SparseMatrix<T, Eigen::ColMajor> SparseMatrixColMajor;
    typedef Eigen::Triplet<Scalar> Triplet;
};

typedef NumTraits<float,Eigen::Dynamic> NumTraitsXf;
typedef NumTraits<double,Eigen::Dynamic> NumTraitsXd;

template <int Top, int Dim>
struct DimTraits {
    enum {top = Top
        dim = Dim};
    typedef DimTraits<top,dim-1> LowerTraits;
    typedef DimTraits<top,dim+1> UpperTraits;
    typedef DimTraits<top,top> TopTraits;
    typedef DimTraits<top,0> BottomTraits;
};
}
#endif

#ifdef _SIMPLICIAL_COMPLEX_H_
#ifndef _SIMPLICIAL_COMPLEX_TRAITS_H_
#define _SIMPLICIAL_COMPLEX_TRAITS_H_
template <typename NT, int D>
class SimplicialComplex;

template <typename NT, typename DT>
class Simplex;
namespace mtao_internal{
template <typename NT, typename DT>
class SimplicialComplexPrivateBase;

template <typename NT, typename DT>
class SimplicialComplexPrivate;


template <typename NT, typename DT>
struct SimplicialComplexTraits {
    typedef NT NumTraits;
    typedef DT DimTraits;
    template <int M=DimTraits::Top>
    struct internal_complex{
        typedef DimTraits<DimTraits::Top,M> dim_traits;
        typedef typename std::conditional<M<=0
        , SimplicialComplexPrivateBase<NT,dim_traits>
        , SimplicialComplexPrivate<NT,dim_traits>
        >::type type;
    };
    template <int M=DimTraits::Top>
    struct simplex{
        typedef DimTraits<DimTraits::Top,M> dim_traits;
        typedef Simplex<NumTraits, dim_traits> type;
    };
    typedef SimplicialComplex<NumTraits,DimTraits::Top> type;
};
};
#endif
#endif


#ifdef _DEC_H_
#ifndef _FORM_TRAITS_H_
#define _FORM_TRAITS_H_

enum class FormType : char {Neither=0, Primal=1,Dual=2, Both=3};
constexpr FormType operator&(FormType a, FormType b) {
    return static_cast<FormType>(static_cast<char>(a) & static_cast<char>(b));
}
namespace mtao_internal{
template <int Dim, typename VectorType, FormType type, int N>
class Form;
template <int Dim_, FormType TypeIn_, int NIn_, FormType TypeOut_, int NOut_, bool isVector_=false>
struct form_operator_traits{
    const static FormType TypeIn = TypeIn_;//Typein = -1 means that this should resolve to a vector
    const static FormType TypeOut = TypeOut_;
    enum {
        Dim = Dim_,
        NIn = NIn_,
        NOut = NOut_,
        isVector = isVector_
    };

};

};

#endif
#ifndef _DEC_TRAITS_H_
#define _DEC_TRAITS_H_


template <typename Complex_, bool Interior_>
class DEC;

namespace mtao_internal{
template <typename Traits, typename Expression>
struct FormExpression;
template <typename Traits, typename Matrixtype>
class FormOperator;
template <int Dim, typename DynamicVector,FormType Type1, int N1>
class Form;

template <typename DECTraits, typename DimTraits>
class OperatorContainerPrivate;
template <typename DECTraits>
class OperatorContainerPrivateBase;
template <typename DECTraits>
class OperatorContainer;

template <typename DECTraits>
class FormFactory;
template <typename Complex_, bool Interior_ = false>
struct dec_traits{
private:
    typedef dec_traits<Complex_, Interior_> Myself;
public:
    typedef Complex_ Complex;
    const static bool Interior = Interior_;
    typedef typename Complex::NumTraits NumTraits;
    typedef typename Complex::DimTraits DimTraits;
    const static int TopDim = DimTraits::Top;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    template <int M=0>
    struct operator_private{
        typedef DimTraits<DimTraits::Top,M> dim_traits;
        typedef typename std::conditional< M>=TopDim
        , OperatorContainerPrivateBase<Myself>
        , OperatorContainerPrivate<Myself,dim_traits>
        >::type type;

        typedef FormOperator<form_operator_traits<TopDim,FormType::Primal,M,FormType::Primal,M+1>,SparseMatrixColMajor> d_primal_type;
        typedef FormOperator<form_operator_traits<TopDim,FormType::Dual,M,FormType::Dual,M+1>, const SparseMatrixColMajor & > d_dual_type;

        typedef FormOperator<form_operator_traits<TopDim,FormType::Dual,M,FormType::Dual,M+1>, const SparseMatrixColMajor > d_dual_interior_type;


        typedef FormOperator<form_operator_traits<TopDim,FormType::Primal,M,FormType::Dual,TopDim-M>, DiagonalMatrix> hodge_primal_type;
        typedef FormOperator<form_operator_traits<TopDim,FormType::Dual,M,FormType::Primal,TopDim-M>, DiagonalMatrix> hodge_dual_type;

        template <FormType Type, bool Interior>
        struct d_type{
            typedef typename std::conditional<Type==FormType::Primal
            , d_primal_type
            , typename std::conditional<Interior, d_dual_interior_type, d_dual_type >::type
            >::type type;
        };
        template <FormType Type>
        struct hodge_type{
            typedef typename std::conditional<Type==FormType::Primal, hodge_primal_type, hodge_dual_type >::type type;
        };

    };
    template <int M=TopDim, FormType Form=FormType::Primal>
    struct dec_operator {
        typedef operator_private<M> operator_type;
        typedef typename operator_private<M>::type type;
        typedef typename std::conditional<(Form==FormType::Primal), typename operator_type::d_primal_type, typename operator_type::d_dual_type>::type d_type;
        typedef typename std::conditional<(Form==FormType::Primal), typename operator_type::hodge_primal_type, typename operator_type::hodge_dual_type>::type h_type;
    };
    template <FormType Type=FormType::Primal, int M=0>
    struct form{
        typedef DimTraits<DimTraits::Top,M> dim_traits;
        typedef Form<DimTraits::Top, typename NumTraits::DynamicVector,Type,M> type;
    };
    typedef OperatorContainer<Myself > operator_type;
    typedef FormFactory<Myself > factory_type;
    typedef DEC<Complex, Interior> type;

};
};
#endif
#endif
