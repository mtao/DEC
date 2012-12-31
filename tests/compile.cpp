#include <iostream>
#include "simplicialComplex.hpp"
#include "mesh.hpp"
int main()
{
    typedef Eigen::Vector3d V;
    std::vector<V> vecs;
    std::vector<mtao::IndexSet<3> > tris;
    for(int i=0; i < 10; ++i)
    {
        vecs.push_back(V(i,i,i));
    }
    tris.push_back(mtao::IndexSet<3>({0,3,5}));
    tris.push_back(mtao::IndexSet<3>({3,0,6}));
    /*
    for(unsigned int i=0; i < 8; ++i)
    {
        tris.push_back(mtao::IndexSet<3>({i,i+1,i+2}));
    }*/

    TriangleMesh sc(tris,vecs);
    auto && verts = sc.Vertices();

    for(auto & m: sc.Vertices())
    {
        std::cout << m.transpose()<< " | ";
    }
   std::cout << std::endl;

    std::cout << sc.Simplices<0>().size() << std::endl;
    for(auto & m: sc.Simplices<0>())
    {
        std::cout << m;
    }
    std::cout << std::endl;
    std::cout << sc.Simplices<1>().size() << std::endl;
    for(auto & m: sc.Simplices<1>())
    {
        std::cout << m;
    }
    std::cout << std::endl;
    std::cout << sc.Simplices<2>().size() << std::endl;
    for(auto & m: sc.Simplices<2>())
    {
        std::cout << m;
    }
    std::cout << std::endl;
    //std::cout << sc.Simplices<3>().size() << std::endl;

    std::cout << verts.size() << std::endl;
    //This file is not included!
    TriangleMesh * mesh = readOBJtoSimplicialComplex<double>("a.obj");
    //std::cout << reinterpret_cast<long>(mesh) << std::endl;
    //std::cout << mesh->Simplices<2>().size() << std::endl;
    return 0;
}
