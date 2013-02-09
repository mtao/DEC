#ifndef DEC_RENDER_H
#define DEC_RENDER_H
#include "dec.hpp"
#include <vector>
namespace mtao{
template <int M, typename SimplicialComplex>
std::vector<unsigned int> simplicesToRenderable(const SimplicialComplex & sc)
{
    static_assert(M > 0 && M <= SimplicialComplex::Dim, "Inappropriate dims asserted");
    auto&& simplices = sc.template constSimplices<M>();
    std::vector<unsigned int> ret(simplices.size() * (M+1));
    for(auto&& s: simplices)
    {
        std::copy(s.getIndexSet().cbegin(),s.getIndexSet().cend(), ret.begin()+s.Index() * (M+1));
        if(s.isNegative())
        {
            unsigned int tmp =  ret[s.Index() * (M+1)+1];
            ret[s.Index() * (M+1)+1] =  ret[s.Index() * (M+1)];
            ret[s.Index() * (M+1)] =  tmp;
        }
    }
    return ret;
}
template <typename Form>
auto formToRenderable(const Form & form) -> std::vector<std::array<float,Form::Traits::NOut+1> > {
    std::vector<std::array<float,Form::Traits::NOut+1> > ret(form.expr.rows());
    for(int i=0; i < form.expr.rows(); ++i) {
        ret[i].fill(form.expr(i));
    }
    return ret;
}
};
#endif
