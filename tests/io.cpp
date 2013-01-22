#include <iostream>
#include "simplicialComplex.hpp"
#include "io.hpp"
#include "dec.hpp"
int main()
{
    TriangleMesh * mesh = readOBJtoSimplicialComplex<double>("a.obj");
    writeOBJfromSimplicialComplex(*mesh,"out.obj");
    return 0;
}
