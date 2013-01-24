#ifndef DEC_OVERLOADS_H
#define DEC_OVERLOADS_H
#ifndef DEC_H
#include "dec.hpp"
#endif

    template <FormType TypeIn, int NIn,FormType TypeOut, int NOut,typename Expr1, typename Expr2>
auto operator+(const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr1> & a, const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr2> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr1>()+std::declval<Expr2>()
            )>
{
    return
        FormExpression<TypeIn, NIn, TypeOut, NOut,
        decltype(
                std::declval<Expr1>()+std::declval<Expr2>()
                )> (a.expr+b.expr);
}

    template <FormType TypeIn, int NIn,FormType TypeOut, int NOut,typename Expr1, typename Expr2>
auto operator-(const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr1> & a, const FormExpression<TypeIn,NIn,TypeOut,NOut,Expr2> & b)
    ->
    FormExpression<TypeIn, NIn, TypeOut, NOut,
    decltype(
            std::declval<Expr1>()-std::declval<Expr2>()
            )>
{
    return
        FormExpression<TypeIn, NIn, TypeOut, NOut,
        decltype(
                std::declval<Expr1>()-std::declval<Expr2>()
                )> (a.expr+b.expr);
}

    template <FormType TypeIn, int NIn,FormType TypeMid,int NMid, FormType TypeOut, int NOut,typename Expr1, typename Expr2>
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

#endif
