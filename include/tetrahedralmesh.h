#ifndef _TETRAHEDRALMESH_H_
#define _TETRAHEDRALMESH_H_

#include "dec.hpp"
/*! \brief Tetrahedral mesh with normals*/

class TetrahedralMesh: public TetrahedralMeshf, public DEC<TetrahedralMeshf>
{
public:
    typedef TetrahedralMeshf SCParent;
    typedef DEC<TetrahedralMeshf> DECParent;
    typedef SCParent::NumTraits NumTraits;
    typedef SCParent::DimTraits DimTraits;
    using SCParent::Dim;
    typedef SCParent::Vector Vector;
    TetrahedralMesh()
        : DECParent(*dynamic_cast<SCParent*>(this))
    {}
    TetrahedralMesh(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
        : SCParent(tuples,vertices)
        , DECParent(*dynamic_cast<SCParent*>(this))
    {
        init();
    }
    TetrahedralMesh(const SCParent & parent)
        : SCParent(parent)
        , DECParent(*dynamic_cast<SCParent*>(this))
    {
        init();
    }

    auto cells() ->decltype(SCParent::simplices<3>()) {return simplices<3>();}
    auto faces() ->decltype(SCParent::simplices<2>()) {return simplices<2>();}
    auto edges() -> decltype(SCParent::simplices<1>())  {return simplices<1>();}
    auto cell(int i) -> decltype(SCParent::simplex<3>(0)) {return simplex<3>(i);}
    auto face(int i) -> decltype(SCParent::simplex<2>(0)) {return simplex<2>(i);}
    auto edge(int i) -> decltype(SCParent::simplex<1>(0)) {return simplex<1>(i);}

    Vector getGradient(int triidx, int tetidx) const;

    ~TetrahedralMesh() {
    }
protected:

    void init();



};
#endif
