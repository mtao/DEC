#include "../../include/trianglemesh.h"
#include "../../include/io.hpp"
#include <vector>


typedef NormalTriangleMesh MeshType;
typedef NormalTriangleMesh::SparseMatrix SparseMatrix;
typedef DEC<MeshType, true > DECType;



template <typename ComplexType>
std::pair<std::vector<typename ComplexType::SCTraits::template simplex<1>::type> , std::vector<typename ComplexType::Vector> > getVoronoiMesh(const std::string& filename) {
    enum {Dim = ComplexType::Dim};
    typedef typename ComplexType::Vector Vector;
    auto mesh = readOBJtoSimplicialComplex<ComplexType>(filename);
    const SparseMatrix & d1 = mesh->template d<1>().expr;

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

    auto&& pr = getVoronoiMesh<MeshType>(argv[1]);

    auto&& vertices = pr.second;
    auto&& edges = pr.first;

    for(auto&& v: vertices) {
        std::cout << "v " << v.transpose() << std::endl;
    }
    for(auto&& e: edges) {
        std::cout << "e " << e[0] << " " << e[1] <<std::endl;
    }

}
