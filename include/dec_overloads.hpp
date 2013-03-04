#ifndef DEC_OVERLOADS_H
#define DEC_OVERLOADS_H
#ifndef DEC_H
#include "dec.hpp"
#endif
#include <typeinfo>

//Notes on deduction
//BOTH_TRAITS is fairly limited in that if an operator has a different storage depending on the trait
//the system will be unable to resolve it
namespace mtao_internal{
template <typename Traits1, typename Traits2, bool isVector>
struct deduceFormTraits
{
    //Traits1 + Traits2
    typedef form_operator_traits<
    (Traits1::Dim == Traits2::Dim) ? Traits1::Dim : -2,
    static_cast<FormType>(Traits1::TypeIn & Traits2::TypeIn),
    (Traits1::NIn == Traits2::NIn) ? Traits1::NIn : -2,
    static_cast<FormType>(Traits1::TypeIn & Traits2::TypeIn),
    (Traits1::NOut == Traits2::NOut) ? Traits1::NOut : -2,
    isVector
    >
    addition;

    //Traits1 * Traits2
    typedef form_operator_traits<
    (Traits1::Dim == Traits2::Dim) ? Traits1::Dim : -2,
    (Traits1::TypeIn & Traits2::TypeOut) ? Traits2::TypeIn : NO_FORM,
    (Traits1::NIn == Traits2::NOut) ? Traits2::NIn : -2,
    (Traits1::TypeIn & Traits2::TypeOut) ? Traits1::TypeOut: NO_FORM,
    (Traits1::NIn == Traits2::NOut) ? Traits1::NOut : -2,
    isVector
    >
    composition;
};
//If either trait type is empty, form is not valid
template <typename Traits>
constexpr bool is_valid_form_trait()
{
    return
            (Traits::Dim >= 0)
            && (Traits::TypeIn | Traits::TypeOut)
            && ((Traits::isVector && Traits::NIn == -1) || Traits::NIn >= 0)
            && (Traits::NOut >= 0)
            ;
}
template <typename Traits1, typename Traits2>
struct deduceFormAdditionTraits
{
private: typedef typename std::enable_if<
    Traits1::isVector == Traits2::isVector,
    typename deduceFormTraits<Traits1,Traits2,Traits1::isVector>::addition >::type Output;
public: typedef typename std::enable_if<is_valid_form_trait<Output>(), Output>::type type;
};
template <typename Traits1, typename Traits2>
struct deduceFormCompositionTraits
{
private: typedef typename deduceFormTraits<Traits1,Traits2,Traits2::isVector>::composition Output;
public: typedef typename std::enable_if<is_valid_form_trait<Output>(), Output>::type type;
};

template <typename Traits1, typename Traits2,typename Expr1, typename Expr2>
auto operator+(const mtao_internal::FormExpression<Traits1,Expr1> & a, const mtao_internal::FormExpression<Traits2,Expr2> & b)
->
const mtao_internal::FormExpression<typename mtao_internal::deduceFormAdditionTraits<Traits1,Traits2>::type,
decltype(
        std::declval<Expr1>()+std::declval<Expr2>()
        )>
{ 
    return
            mtao_internal::FormExpression<typename mtao_internal::deduceFormAdditionTraits<Traits1,Traits2>::type,
            decltype(
                std::declval<Expr1>()+std::declval<Expr2>()
                )> {a.expr+b.expr};
}
template <typename Traits1, typename Traits2,typename Expr1, typename Expr2>
auto operator-(const mtao_internal::FormExpression<Traits1,Expr1> & a, const mtao_internal::FormExpression<Traits2,Expr2> & b)
->
const mtao_internal::FormExpression<typename mtao_internal::deduceFormAdditionTraits<Traits1,Traits2>::type,
decltype(
        std::declval<Expr1>()-std::declval<Expr2>()
        )>
{
    return
            mtao_internal::FormExpression<typename mtao_internal::deduceFormAdditionTraits<Traits1,Traits2>::type,
            decltype(
                std::declval<Expr1>()-std::declval<Expr2>()
                )> {a.expr-b.expr};
}
template <typename Traits1, typename Traits2,typename Expr1, typename Expr2>
auto operator*(const mtao_internal::FormExpression<Traits1,Expr1> & a, const mtao_internal::FormExpression<Traits2,Expr2> & b)
->
const mtao_internal::FormExpression<typename mtao_internal::deduceFormCompositionTraits<Traits1,Traits2>::type,
decltype(
        std::declval<Expr1>()*std::declval<Expr2>()
        )>
{ 
    return
            mtao_internal::FormExpression<typename mtao_internal::deduceFormCompositionTraits<Traits1,Traits2>::type,
            decltype(
                std::declval<Expr1>()*std::declval<Expr2>()
                )> {a.expr*b.expr};
}
constexpr int toPrimalDim(FormType type, int N, int TopD){return 0;}
};
/*
    template <
        FormType TypeIn1, int NIn1,FormType TypeOut1, int NOut1, typename Expr1,
        FormType TypeIn1, int NIn1,FormType TypeOut1, int NOut1,typename Expr2
        >
auto operator-(const FormExpression<TypeIn1,NIn1,TypeOut1,NOut1,Expr1> & a, const FormExpression<TypeIn2,NIn2,TypeOut2,NOut2,Expr2> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr1>()-std::declval<Expr2>(:)
            )>
{
    static_assert(
                ()
                );
    return
        FormExpression<TypeIn, NIn, TypeOut, NOut,
        decltype(
                std::declval<Expr1>()-std::declval<Expr2>()
                )> (a.expr-b.expr);
}

    template <FormType TypeIn1, int NIn1,FormType TypeOut1,int NOut1, FormType TypeOut2, int NOut2,typename Expr2>
auto operator*(const FormExpression<TypeMid,NMid,TypeOut,NOut,Expr1> & a, const FormExpression<TypeIn,NIn,TypeMid,NMid,Expr2> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr1>()+std::declval<Expr2>()
            )>
{
    return FormExpression<TypeIn, NIn, TypeOut, NOut,
           decltype(
                   std::declval<Expr1>()+std::declval<Expr2>()
                   )> (a.expr*b.expr);
}

    template <typename Scalar,FormType TypeIn, int NIn,FormType TypeOut, int NOut,typename Expr>
auto operator*(typename std::enable_if<std::is_scalar<Scalar>::value,Scalar>::type a, const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Scalar>()+std::declval<Expr>()
            )>
{
    return FormExpression<TypeIn, NIn, TypeOut, NOut,
           decltype(
                   std::declval<Scalar>()+std::declval<Expr>()
                   )> (a*b.expr);
}
    template <typename Scalar, FormType TypeIn, int NIn,FormType TypeOut, int NOut,typename Expr>
auto operator*(const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr> & a, typename std::enable_if<std::is_scalar<Scalar>::value,Scalar>::type b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr>()+std::declval<Scalar>()
            )>
{
    return FormExpression<TypeIn, NIn, TypeOut, NOut,
           decltype(
                   std::declval<Expr>()+std::declval<Scalar>()
                   )> (a*b.expr);
}
*/
#endif
