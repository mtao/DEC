#ifndef _SIMPLEX_LIBRARY_H_
#define _SIMPLEX_LIBRARY_H_
#include "library.hpp"

template <typename NT, typename DT>
struct Simplex;

namespace DEC {
    namespace internal {
        template <typename NT, typename DT, int N=DT::dim>
            struct Library<Simplex<NT,DT> >::object<N> {
                typedef Simplex<NT,DimTraits<DT::top,N> > type;
            };
    }
}

#endif
