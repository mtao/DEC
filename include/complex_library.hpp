#ifndef _COMPLEX_LIBRARY_H_
#define _COMPLEX_LIBRARY_H_
#include "library.hpp"


template <typename ObjectType>
Complex;

namespace DEC{ namespace internal{ 
    template <typename NumTraits>
        class ComplexBase_;
    template <typename ObjectType>
        class Complex_;

}}

namespace DEC{namespace internal{
    template <typename ObjectType>
        struct Library<Complex<ObjectType> > {
            typedef Complex<ObjectType> object_type;
            typedef ObjectType::num_traits num_traits;
            typedef ObjectType::dim_traits dim_traits;
            template <int M=dim_traits::top>
                struct internal_complex{
                    typedef DimTraits<dim_traits::top,M> dim_traits;
                    typedef typename std::conditional<M<=0
                        , ComplexBase_<num_traits>
                        , Complex_<ObjectType::object<M>::type >
                        >::type type;
                };
            template <int M=dim_traits::top>
                struct internal_object{
                    typedef DimTraits<dim_traits::top,M> dim_traits;
                    typedef Library<ObjectType>::object<M> type;
                };
                template <int N=dim_traits::dim>
                    struct object {
                        typedef Complex<internal_object<N>::type> type;
                    };
        };
};
}}
#endif

