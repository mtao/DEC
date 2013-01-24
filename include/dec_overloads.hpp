#ifndef DEC_OVERLOADS_H
#define DEC_OVERLOADS_H
#ifndef DEC_H
#include "dec.hpp"
#endif

template <FormType Type1,int N1,typename Expr1, FormType Type2, int N2, typename Expr2>
auto operator+(const FormExpression<Type1,N1,Expr1> & a, const FormExpression<Type2,N2,Expr2> & b)
->
FormExpression<Type1, N1,
decltype(
        std::declval<Expr1>()+std::declval<Expr2>()
        )>
{
    static_assert(Type1 == Type2, "Addition can't combine different forms");
    static_assert(N1 == N2, "Addition can't combine different forms");
    return
            (static_cast<Expr1>(a)+static_cast<Expr2>(b));
}



#endif
