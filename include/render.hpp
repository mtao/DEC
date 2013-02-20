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
constexpr unsigned int MeshSize() {
    return (Form::Traits::TypeOut == PRIMAL_FORM) *
        Form::Traits::NOut+1 +
    (Form::Traits::TypeOut == DUAL_FORM) *
        Form::Traits::Dim-Form::Traits::NOut+1;
}
template <typename Form>
constexpr bool shouldFormGenUseIndices() {
    return (Form::Traits::NOut == 2 && Form::Traits::TypeOut == DUAL_FORM);
}

template <typename Form>
auto formToRenderable(const Form & form
                      , const std::vector<unsigned int> & indices = std::vector<unsigned int>()
        ) -> std::vector<std::array<float,MeshSize<Form>() > > {
    std::vector<
            std::array<float,MeshSize<Form>() >
            > ret(shouldFormGenUseIndices<Form>()
                                                    ?
                                                      indices.back()
                                                    :
                                                        form.expr.rows());
    if(!shouldFormGenUseIndices<Form>()) {
    for(int i=0; i < form.expr.rows(); ++i) {
        ret[i].fill(form.expr(i));
    }
    } else {
        unsigned int j=0;
        for(int i=0; i < indices.size(); ++i) {
            const unsigned int maxind = indices[i];
            for(;j<maxind; ++j) {
                ret[j].fill(form.expr(i));
            }
        }
    }
    return ret;
}
    /*
template <typename SimplicialComplex, typename Form>
auto formToFluxRenderable(const SimplicialComplex & sc, const Form & form) -> std::vector<std::array<float,Form::Traits::NOut+1> > {

    static const int N = SimplicialComplex::Dim;
    auto&& nsimplices = sc.simplices();
    auto&& n1simplices = sc.template simplices<N-1>();
    std::vector<std::array<float,Form::Traits::NOut+1> > ret(nsimplices.size());
    std::transform(ind.cbegin(), ind.cend(), nsimplices.begin(), [&form](const Indices & ind) {
        for(int i=0; i <=N; ++i) {
        }

    });
    for(int i=0; i < form.expr.rows(); ++i) {
        ret[i].fill(form.expr(i));
    }
    return ret;
}
    */
};
#endif
