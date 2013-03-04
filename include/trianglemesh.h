#ifndef _TRIANGLEMESH_H_
#define _TRIANGLEMESH_H_

#include "dec.hpp"
/*! \brief Triangle mesh with normals*/

class NormalTriangleMesh: public TriangleMeshf, public DEC<TriangleMeshf>
{
public:
    enum NormalType {Equal_Normal, Area_Normal, Angle_Normal, Mean_Curvature_Normal, Sphere_Normal};
    enum AdvectionType {Semilagrangian_Advection};
    typedef TriangleMeshf SCParent;
    typedef DEC<TriangleMeshf> DECParent;
    typedef SCParent::NumTraits NumTraits;
    typedef SCParent::DimTraits DimTraits;
    typedef typename DECParent::Nm1Form VelocityFormType;
    using SCParent::Dim;
    typedef SCParent::Vector Vector;
    typedef SCParent::Scalar Scalar;
    typedef SCParent::SparseMatrixColMajor SparseMatrixColMajor;
    NormalTriangleMesh()
        : DECParent(*dynamic_cast<SCParent*>(this))
    {}
    NormalTriangleMesh(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
        : SCParent(tuples,vertices)
        , DECParent(*dynamic_cast<SCParent*>(this))
    {
        init();
    }
    NormalTriangleMesh(const SCParent & parent)
        : SCParent(parent)
        , DECParent(*dynamic_cast<SCParent*>(this))
    {
        init();
    }

    auto faces() ->decltype(SCParent::simplices<2>()) {return simplices<2>();}
    auto edges() -> decltype(SCParent::simplices<1>())  {return simplices<1>();}
    auto face(int i) -> decltype(SCParent::simplex<2>(0)) {return simplex<2>(i);}
    auto edge(int i) -> decltype(SCParent::simplex<1>(0)) {return simplex<1>(i);}

    ~NormalTriangleMesh() {
    }
    std::vector<Vector> getNormals(NormalType type = Equal_Normal);
    void advect(VelocityFormType & form, Scalar dt, AdvectionType type = Semilagrangian_Advection);
protected:
    std::vector<Vector> m_face_normals;

    std::vector<Vector> equalNormal();
    std::vector<Vector> areaNormal();
    std::vector<Vector> angleNormal();
    std::vector<Vector> meanCurvatureNormal();
    std::vector<Vector> sphereNormal();
    void semilagrangianAdvection(VelocityFormType & form, Scalar dt);
    void init();
    void computeNormals();



};
#endif
