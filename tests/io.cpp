#include <iostream>
#include "simplicialComplex.hpp"
#include "io.hpp"
#include "dec.hpp"
int main()
{
    typedef Eigen::Vector3d V;
    std::vector<V> vecs;
    std::vector<mtao::IndexSet<3> > tris;
    for(int i=0; i < 10; ++i)
    {
        vecs.push_back(V(i,i,i));
    }
    vecs[0] = V(0,0,0);
    vecs[3] = V(1,1,0);
    vecs[5] = V(1,0,0);
    vecs[6] = V(0,1,0);
    tris.push_back(mtao::IndexSet<3>({0,3,5}));
    tris.push_back(mtao::IndexSet<3>({3,0,6}));

    TriangleMesh sc(tris,vecs);
    writeOBJfromSimplicialComplex(sc,"out.obj");
    return !readOBJtoSimplicialComplex<TriangleMesh>("out.obj");
}
