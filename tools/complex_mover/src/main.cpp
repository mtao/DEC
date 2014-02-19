#include <iostream>
#include <fstream>
#include <openvdb/openvdb.h>
#include <openvdb/tools/GridOperators.h>
#include <openvdb/tools/LevelSetAdvect.h>
#include <openvdb/tools/ValueTransformer.h>
#include <openvdb/tools/MeshToVolume.h>
#include <sstream>
#include <iomanip>
#include "trianglemesh.h"
#include "io.hpp"
#include <vector>
#include "../include/conversion.h"
#include <openvdb/math/Stencils.h>


typedef NormalTriangleMesh MeshType;
typedef NormalTriangleMesh::SparseMatrix SparseMatrix;
typedef DEC<MeshType, true > DECType;








template <typename ComplexType>
std::pair<std::vector<typename ComplexType::SCTraits::template simplex<1>::type> , std::vector<typename ComplexType::Vector> > getVoronoiMesh(const typename ComplexType::Ptr& mesh) {
    enum {Dim = ComplexType::Dim};
    typedef typename ComplexType::Vector Vector;
    const SparseMatrix & d1 = mesh->template d<Dim-1>().expr;

    typedef typename ComplexType::SCTraits::template simplex<1>::type EdgeType;
    std::vector<EdgeType> indices(mesh->template simplices<Dim-1>().size());
    std::vector<Vector> vertices(mesh->template simplices<Dim>().size());


    auto&& faces = mesh->template simplices<Dim>();
    typedef typename ComplexType::NSimplex FaceType;

    std::transform(faces.begin(), faces.end(), vertices.begin(), [&](const FaceType& face) -> Vector {
            return face.Center();
            });


    for(int i=0; i < mesh->template simplices<Dim-1>().size(); ++i) {
        auto&& s = mesh->template simplex<Dim-1>(i);
        SparseMatrix::InnerIterator it(d1, s.Index());
        auto&& e = indices[i];
        if(!it) continue;
        e[0] = it.row();
        ++it;
        if(!it) {
            auto&& v = s.Center();
            vertices.push_back(v);
            e[1] = vertices.size()-1;
        } else {
            e[1] = it.row();
        }

    }
    return std::make_pair(indices,vertices);

}

int main(int argc, char * argv[]) {
    //auto mesh = readOBJtoSimplicialComplex<MeshType>(argv[1]);
    openvdb::initialize();

    MeshType::Ptr mesh(readOBJtoSimplicialComplex<MeshType>(argv[1]));
    auto&& sc = *mesh;
    std::vector<openvdb::math::Vec3s> verts(sc.vertices().size());
    std::vector<openvdb::Vec4I> poly(sc.simplices().size());

    typedef MeshType::Vector Vector;
    typedef MeshType::NSimplex Simplex;

    std::transform(sc.vertices().begin(),sc.vertices().end(),verts.begin(),[&](const Vector& v) -> openvdb::math::Vec3s {
            return openvdb::math::Vec3s(v[0],v[1],v[2]);
            });
    std::transform(sc.simplices().begin(),sc.simplices().end(),poly.begin(), [&](const Simplex& s) -> openvdb::Vec4I {
            if(s.isNegative()) {
            return openvdb::Vec4I(s[0],s[2],s[1],openvdb::util::INVALID_IDX);
            } else {
            return openvdb::Vec4I(s[0],s[1],s[2],openvdb::util::INVALID_IDX);
            }
            });
    float scale = 1.0 / 10;


    auto&& tptr = openvdb::math::Transform::createLinearTransform(scale);
    std::cout << "Creating levelset" << std::endl;
    auto&& gptr = openvdb::tools::meshToLevelSet<openvdb::FloatGrid>(*tptr, verts, poly);
    auto&& grid = *gptr;
    typedef openvdb::FloatGrid GridType;

    typedef typename openvdb::tools::Gradient<GridType> Gradient;
    Gradient grad(grid);
    typedef typename Gradient::OutGridType GradGrid;
    typename GradGrid::Ptr gradGrid = grad.process();
    typedef typename openvdb::tools::DiscreteField<GradGrid> GradDiscreteField;
    GradDiscreteField discField(*gradGrid);


    auto&& pr = getVoronoiMesh<MeshType>(mesh);

    auto&& dual_vertices = pr.second;
    auto&& dual_edges = pr.first;
    auto&& vertices = sc.vertices();
    auto&& edges = sc.template simplices<1>();

    typedef openvdb::math::BoxStencil<GridType> BoxStencil;
    BoxStencil bs(grid);

    for(auto&& v: vertices) {
        auto&& vdbv = eigen2vdb(v);
        auto&& iv = grid.worldToIndex(vdbv);
        bs.moveTo(iv);
        auto&& vec = bs.interpolate(iv);
        v = v - vdb2eigen(vec);

    }
    for(auto&& v: vertices) {
        std::cout << "v " << v.transpose() << std::endl;
    }
    for(auto&& v: dual_vertices) {
        std::cout << "v " << v.transpose() << std::endl;
    }
    for(auto&& e: edges) {
        std::cout << "e " << e[0] << " " << e[1] << std::endl;
    }
    for(auto&& e: dual_edges) {
        std::cout << "e " << e[0]+vertices.size() << " " << e[1]+vertices.size() << std::endl;
    }

}

