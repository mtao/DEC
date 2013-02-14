#ifndef DEC_KDTREE_H
#define DEC_KDTREE_H


#include <memory>
#include <Eigen/Geometry>
#include "simplex.hpp"
template <typename Vector>
class KDTree;
namespace mtao {

template <typename Object>
struct KDNode{
    //const KDTree<Vector> & tree;
    typedef typename Object::Vector Vector;
    Eigen::AlignedBox<Vector::Scalar, Vector::Dim> bbox;
    int index = -1;
    char dim = -1;
    std::unique_ptr<KDNode<Object> > left;
    std::unique_ptr<KDNode<Object> > right;
    KDNode(const KDTree<Vector> & tree) :tree(tree) {}
    //this constructor allows for implit conversion from vector -> kdnode which makes kdnode / vector comparison possible
    KDNode(const Vector & v, int dim_=-1, int index_=-1): bbox(v,v), dim(dim_), index(index_) {}
    KDNode(const Object & obj, int dim_=-1, int index_=-1): bbox(obj.boundingBox()), dim(dim_), index(index_) {}
};
template <typename Polyhedron, int N>
struct KDPolyhedronBBoxContainer{
    typedef Polyhedron::NumTraits::Vector Vector;
    const decltype(m_bbox) & boundingBox() {
        return m_bbox;
    }
    KDPolyhedronBBoxContainer(const IndexSet<N> & indices, const std::vector<Vector> & verts) {
        for(int i=0; i < N; ++i)
        {
            m_bbox = m_bbox.merged(
                        decltype(m_bbox)(verts[indices[i]],verts[indices[i]])
                    );
        }
    }
    private:
    Eigen::AlignedBox<Vector::Scalar, Vector::Dim> m_bbox;
};

template <typename Vector>
bool operator<(const KDNode<Vector> & a, const KDNode<Vector> & b) {
    int dim=a.dim;
    if(dim < 0) dim = b.dim
    return a.tree.vertices[a.index](dim) < b.tree.vertices[b.index](dim);
}

//assumes a persisting vector of Vectors that the tree represents
//which is pretty sane considering the KDTree only represents a static set of verts...
template <typename Vector>
class KDTree{
public:
    typedef typename Vector::Scalar Scalar;
    typedef KDNode<Vector> Node;
    friend class Node;
    KDTree() {}
    KDTree(const std::vector<Vector> & vertices, bool doprocess = true)
        : vertices(vertices)
        , nodes(vertices.size(),Node(*this)){
        if(doprocess) {
            process();
        }
    }
    int findIndex(const Vector & vector, int idx) const;
    int findIndex(const Vector & vector) const {
        return findIndex(vector,0);
    }
    const Node & findNode(const Vector & vector) const {
        return nodes[findIndex(vector)];
    }


protected:
    void process();//mainly here so simplexkdtree can populate the simplices in time...
    void insert(int idx);
    std::vector<Vector> & vertices;
    std::vector<Node> nodes;
private:
    int lastUnfilledPosition = 0;
};

/*
template <typename SimplicialComplex_>
class SimplexKDTree: public KDTree<typename SimplicialComplex_::Traits::Vector> {
public:
    typedef SimplicialComplex_ SimplicialComplex;
    typedef typename SimplicialComplex::Traits::Vector Vector;
    typedef typename SimplicialComplex::Traits::Scalar Scalar;
    typedef KDTree<Vector> Parent;

    SimplexKDTree(const SimplicialComplex & sc)
        : Parent(centers,false)
        , centers(sc.numSimplices())
        , radii(sc.numSimplices()){
        auto&& simplices = sc.simplices();
        for(int i=0; i < centers.size(); ++i) {
            centers[i] = simplices[i].Center();
            radii[i] =
                    (sc.vertices()[simplices[i][0]] - centers[i]).norm();
        }
        Parent::process();


    }

private:
    std::vector<Vector> centers;
    std::vector<Scalar> radii;
};
*/


template <typename Vector>
void KDTree::Vector>::process()
{
    nodes.resize(vertices.size(), Node(*this));
    for(int i=0; i < vertices.size(); ++i) {
        insert(v);
    }
}
template <typename Vector>
int KDTree<Vector>::findIndex(const Vector & vertex, int idx) {
    while(true) {
        Node& n = nodes[idx];
        if(vertex < n) {
            if(n.left == -1) {
                return idx;
            } else {
                return findIndex(n.left);
            }
        } else if(n.right == -1) {
            return idx;
        } else {
            return findIndex(n.right);
        }


    }
}


template <typename Vector>
void KDTree<Vector>::insert(int idx) {
    auto&& vert = vertices[idx];
    Node & parent = nodes[findIndex(vert)];
    int newind = lastUnfilledPosition++;
    Node& me = nodes[newind];
    me.index = idx;
    me.dim = (parent.dim+1)%vert.rows();
    if(vert < parent) {
        parent.left = newind;
    } else {
        parent.right = newind;
    }
}

#endif
