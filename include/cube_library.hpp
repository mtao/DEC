#ifndef _CUBE_LIBRARY_H_
#define _CUBE_LIBRARY_H_
#include "library.hpp"

template <typename NT, typename DT>
struct Cube;

namespace DEC {
    namespace internal {
        template <typename NT, typename DT, int N=DT::dim>
            struct Library<Cube<NT,DT> >::object<N> {
                typedef Cube<NT,DimTraits<DT::top,N> > type;
            };
    }
}

#endif
