#ifndef DEC_H
#define DEC_H
#include "simplicialComplex.hpp"

enum FormType {PRIMAL, DUAL};
template <typename SimplicialComplex,FormType Type, int N>
class Form: public SimplicialComplex::NumTraits::DynamicVector
{
public:
    typedef typename SimplicialComplex::NumTraits NumTraits;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::DynamicVector Parent;
    Form(): Parent(0) {}
    Form(const SimplicialComplex & sc)
        : Parent(
              (Type == PRIMAL) ?  sc.template numSimplices<N>() : sc.template numSimplices <SimplicialComplex::Dim-N>()
                                  )
    {
        static_assert(N <= SimplicialComplex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
        this->setConstant(Scalar(0));
    }
    void init(const SimplicialComplex & sc)
    {
        resize((Type == PRIMAL) ?  sc.template numSimplices<N>() :sc.template numSimplices<SimplicialComplex::Dim-N>()
                                   );
        setConstant(Scalar(0));

    }
};

#endif
