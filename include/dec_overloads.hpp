#ifndef DEC_OVERLOADS_H
#define DEC_OVERLOADS_H
#ifndef DEC_H
#include "dec.hpp"
#endif
#include <typeinfo>

template <FormType F1, FormType F2>
struct deduceFormType
{
    typedef form_traits<(
            Traits1::TypeIn==BOTH_FORM)?Traits2:()
        type;
};

    template <typename Traits1, typename Traits2,typename Expr1, typename Expr2>
auto operator+(const FormExpression<Traits1,Expr1> & a, const FormExpression<Traits2,Expr2> & b)
    ->
    const FormExpression<deduceTraits_addition<Traits1,Traits2>::type,
    decltype(
            std::declval<Expr1>()+std::declval<Expr2>()
            )>
{
    static_assert(deduceType(Traits1,Traits2), "Expression types are not compatible");
    return
        FormExpression<deduceTraits_addition<Traits1,Traits2>::type,
        decltype(
                std::declval<Expr1>()+std::declval<Expr2>()
                )> {a.expr+b.expr};
}
constexpr int toPrimalDim(FormType type, int N, int TopD){return 0;}
/*
    template <
        FormType TypeIn1, int NIn1,FormType TypeOut1, int NOut1, typename Expr1,
        FormType TypeIn1, int NIn1,FormType TypeOut1, int NOut1,typename Expr2
        >
auto operator-(const FormExpression<TypeIn1,NIn1,TypeOut1,NOut1,Expr1> & a, const FormExpression<TypeIn2,NIn2,TypeOut2,NOut2,Expr2> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr1>()-std::declval<Expr2>()
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
