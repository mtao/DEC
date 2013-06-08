#ifndef _COMPLEX_H_
#define _COMPLEX_H_


//Stores a complex, assumes that the Simplex has all of the necessary information
template <typename ObjectType>
class Complex: public mtao_internal::ComplexPrivate<ObjectType,mtao_internal::DimTraits<N_,N_> > {
private: typedef mtao_internal::ComplexPrivate<NT,mtao_internal::DimTraits<N_,N_> > PrivateParent;
public:
         typedef ObjectType object_type;
    typedef NT num_traits;
    typedef mtao_internal::DimTraits<N_,N_> dim_traits;
    typedef mtao_internal::ComplexTraits<num_traits, dim_traits> SCTraits;
    static const int Dim = N_;
    static const int N = N_;
    static const int TopDim = N_;
    static const int EmbeddedDim = NT::Dim;
    PLAIN_CLASS_NUM_DEFS
    typedef ObjectType::trait_library object_trait_library;

    template <int M=N>
    struct TraitsContainer{
        typedef mtao_internal::DimTraits<dim_traits::Top,M> dim_traits;
        typedef typename std::conditional<M==0
        , mtao_internal::ComplexPrivateBase<NT,typename dim_traits::BottomTraits >
        , mtao_internal::ComplexPrivate<NT,dimtype >
        >::type complextype;
        typedef Simplex<NT,dimtype> simplextype;

    };
    //====================================================
    //====================================================
    //==               Constructors                     ==
    //====================================================
    //====================================================
    Complex() {}
    Complex(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices)
        : PrivateParent(simplices,vertices)
    {
    }
    Complex(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }
    Complex(const std::vector<std::array<unsigned int, N+1> > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }


    Complex(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }

    //====================================================
    //====================================================
    //==               Simplex Accessors                ==
    //====================================================
    //====================================================
    template <int M=Dim>
    size_t numSimplices()const{return TraitsContainer<M>::complextype::m_simplices.size();}

    template <int M=Dim>
    typename TraitsContainer<M>::complextype::NSimplex & simplex(int ind) {
        return TraitsContainer<M>::complextype::simplexByIndex(ind);
    }
    template <int M=Dim>
    const typename TraitsContainer<M>::simplextype & simplex(int ind) const {
        return TraitsContainer<M>::complextype::simplexByIndex(ind);
    }

    template <int M=Dim>
    std::vector< typename TraitsContainer<M>::simplextype > & simplices() {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return TraitsContainer<M>::complextype::m_simplices;
    }
    template <int M=Dim>
    const std::vector<typename TraitsContainer<M>::simplextype> & constSimplices() const {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return TraitsContainer<M>::complextype::m_simplices;
    }


    template <int M=Dim>
    const std::map<int,int> & indexToSimplex() const {return TraitsContainer<M>::complextype::m_indexToSimplex;}

    //====================================================
    //====================================================
    //==              Topology Accessors                ==
    //====================================================
    //====================================================
    template <int M=Dim>
    const SparseMatrix & b() const{return TraitsContainer<M>::complextype::m_boundary;}
    template <int M=Dim>
    const DiagonalMatrix & interior() const{return TraitsContainer<M>::complextype::m_interior;}
    //====================================================
    //====================================================
    //==              Vertex Accessors                  ==
    //====================================================
    //====================================================

    /*! Accessors for the vertices*/
    std::vector<Vector> & vertices(){return TraitsContainer<0>::complextype::m_vertices;}
    const std::vector<Vector> & constVertices()const{return TraitsContainer<0>::complextype::m_vertices;}
    Vector & vertex(unsigned int ind) { return TraitsContainer<0>::complextype::m_vertices[ind];}
    const Vector & vertex(unsigned int ind) const{ return TraitsContainer<0>::complextype::m_vertices[ind];}
    template <int M=Dim>
    auto vertices(const typename TraitsContainer<M>::simplextype & s) const
    -> decltype(TraitsContainer<M>::complextype::verticesBySimplex(s)) {
        return TraitsContainer<M>::complextype::verticesBySimplex(s);
    }


    //====================================================
    //====================================================
    //==            Whitney Accessors                   ==
    //====================================================
    //====================================================
    template <int M=Dim>
    const typename TraitsContainer<M>::complextype::WhitneyBasis & whitneyBasis(unsigned int ind) const {
        return TraitsContainer<M>::complextype::m_whitneyBases[ind];
    }
    template <int M=Dim>
    auto whitneyBasis(const typename TraitsContainer<M>::simplextype & simplex) const
    -> decltype(TraitsContainer<M>::complextype::m_whitneyBases[simplex.Index()]) {
        return TraitsContainer<M>::complextype::m_whitneyBases[simplex.Index()];
    }
    template <int M=Dim>
    auto whitneyCenters(const typename TraitsContainer<M>::simplextype & simplex) const
    -> decltype(TraitsContainer<M>::complextype::m_whitneyCenters[simplex.Index()]) {
        return TraitsContainer<M>::complextype::m_whitneyCenters[simplex.Index()];
    }
    template <int M=Dim>
    auto barycentricCoords(const typename TraitsContainer<M>::simplextype & s, const Vector & v) const
    -> typename TraitsContainer<M>::complextype::BarycentricCoordinates {
        return TraitsContainer<M>::complextype::barycentricCoords(s,v);
    }
    template <int M=Dim>
    bool inSimplex(const NSimplex & s, const Vector & v) const {
        return TraitsContainer<M>::complextype::inSimplex(s,v);
    }
    template <int M=Dim>
    void projectToSimplexInPlace(const NSimplex & s, Vector & v) const {
        TraitsContainer<M>::complextype::projectToSimplexInPlace(s,v);
    }

    template <int M=Dim>
    Vector projectToSimplex(const NSimplex & s, const Vector & v) const {
        Vector nv = v;
        TraitsContainer<M>::complextype::projectToSimplexInPlace(s,nv);
        return nv;
    }

};


typedef Complex<mtao_internal::num_traits<float, 3>, 2> TriangleMeshf;
typedef Complex<mtao_internal::num_traits<float, 3>, 3> TetrahedralMeshf;
typedef Complex<mtao_internal::num_traits<double, 3>, 2> TriangleMeshd;
typedef Complex<mtao_internal::num_traits<double, 3>, 3> TetrahedralMeshd;
//Default to double :)
typedef Complex<mtao_internal::num_traits<double, 3>, 2> TriangleMesh;
typedef Complex<mtao_internal::num_traits<double, 3>, 3> TetrahedralMesh;



