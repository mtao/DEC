#include <iostream>
#include "simplicialComplex.hpp"
#include "io.hpp"
#include "dec.hpp"
#include "render.hpp"
int main()
{
    std::cout << "Running compile test: " << std::endl;
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
    std::cout << "Vertices: " << sc.vertices().size() << std::endl;
    for(auto & m: sc.vertices())
    {
        std::cout << m.transpose()<< " | ";
    }
    std::cout << std::endl;

    std::cout << sc.simplices<0>().size() << std::endl;
    for(auto & m: sc.simplices<0>())
    {
        std::cout << m;
    }
    std::cout << std::endl;

    std::cout << sc.simplices<1>().size() << std::endl;
    for(auto & m: sc.simplices<1>())
    {
        std::cout << m;
    }
    std::cout << std::endl;
    auto simps = mtao::simplicesToRenderable<1>(sc);
    for(auto&& a: simps)
    {
        std::cout << a << " ";
    }
    std::cout << std::endl;
    std::cout << sc.simplices<2>().size() << std::endl;
    for(auto & m: sc.simplices<2>())
    {
        std::cout << m;
    }
    std::cout << std::endl;
    simps = mtao::simplicesToRenderable<2>(sc);
    for(auto&& a: simps)
    {
        std::cout << a << " ";
    }
    std::cout << std::endl;
    writeSimplicialComplextoStream(sc,std::cout);
    std::cout << "Interior calcs: " << std::endl;

    std::cout << sc.template interior<2>().diagonal().transpose() << std::endl;
    std::cout << sc.template interior<1>().diagonal().transpose() << std::endl;
    std::cout << sc.template interior<0>().diagonal().transpose() << std::endl;


    std::cout << "DEC operators" << std::endl;
    DEC<TriangleMesh> dec(sc);
    std::cout << dec.template d<0>().constData() << std::endl;
    std::cout << dec.template d<1>().constData() << std::endl;
    /*
    std::cout << (dec.template d<1>() * dec.template d<0>()) << std::endl;
    std::cout << (dec.template d<1>() * dec.template d<0>()).norm() << std::endl;
    */
    //std::cout << one.transpose().squaredNorm() << std::endl;
    return 0;
}
