DEC
===

My DEC (Discrete Exterior Calculus) implementation which now marginally works.  To some degree this was an experiment in how much I could stretch C++ to allow myself to write code "mathematically".

##Requirements
* gcc4.7 or clang3.1 (or some compiler with decent c++11 support)
* cmake
* Eigen, which should be 3.1 but might require the development branch

## Features
* This implements the core features required in a DEC implementation such as definitions for the discrete exterior derivative and the hodge star.
* Compiles really really slowly  due to all of the type deduction I require the compiler to deal with.
* The simplicial complex (nesh) representation follows an is-a relationship fairly strictly in that simplicial complexes of n dimensions inherit from simplicial complexes of m dimensions, for all m < n.
* Typechecking for differential forms and operators.
* Can compose differential operators using the syntax

    h(d(h(d<1>()))) + d(h(d(h<1>())))

to create an expression template that sits ontop of Eigen's expresison template system.  This is all done through the above mentioned typechecking and a template deduction.
