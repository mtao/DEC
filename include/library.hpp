#ifndef _LIBRARY_H_
#define _LIBRARY_H_
namespace DEC{
    namespace internal {
        template <ObjectType>
            struct Library{
                typedef ObjectType object_type;
                typedef NT num_traits;
                typedef DT dim_traits;
                template <int N=dim_traits::dim>
                    struct object;

        };
    }
}

#endif
