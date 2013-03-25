#ifndef _SIMPLICIAL_COMPLEX_H_
#define _SIMPLICIAL_COMPLEX_H_

#include <set>
#include "simplex.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <iostream>
#include "util.hpp"

//Input N-simplices to generate the whole simplicial complex
//Assumes that the input mesh of simplices is manifold

//NumTraits::N identifies the dimension of the space
//N determines the maximal simplex dimension

namespace mtao_internal
{
/*! A 0-simplicial complex only manages
* the vertices and some trivalish set of simplices.  vertices without any higher order
* information are considered bad data and though stored here, they do not have simplices
* associated with them.
*/
template <typename NT, typename DT>
class SimplicialComplexPrivateBase
{
protected:
    template <int M>
    struct TraitsContainer{
        typedef typename std::enable_if<M==0
        , SimplicialComplexPrivateBase<NT,mtao_internal::dimensional_traits<DT::Top, M> > >::type type;
    };
    typedef NT NumTraits;
    static const int Dim = 0;
    static const int TopDim = DT::Top;
    typedef DT DimTraits;
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef Simplex<NumTraits,DT> NSimplex;

    //====================================================
    //====================================================
    //==               Constructors                     ==
    //====================================================
    //====================================================

    /*! Empty Constructor*/
    SimplicialComplexPrivateBase() {}

    /*! Constructor that just takes vertices and creates no simplices.
      * This shouldn't be called by anything except parent simplicial complex classes which will generate the simplices for us.
      */
    SimplicialComplexPrivateBase(const std::vector<Vector> & vertices)
    {
        m_vertices=vertices;
    }

    /*! This constructor shouldn't ever really be called because there's no point in a 0-dimensional simplicial complex...
      */
    SimplicialComplexPrivateBase(const std::vector<NSimplex >& simplices, const std::vector<Vector> & vertices)
    {
        m_simplices = simplices;
        m_vertices=vertices;
    }






    //====================================================
    //====================================================
    //==            Simplex Ordering Helpers            ==
    //====================================================
    //====================================================
    /*! Provides an external accessor for information about the locations of simplices ordered by index*/
    template <int M=0>
    const std::map<int,int> & indexToSimplex() const {return TraitsContainer<M>::complextype::m_indexToSimplex;}
    template <int M=0>
    int indexToSimplex(int ind) const {return TraitsContainer<M>::complextype::m_indexToSimplex[ind];}
    NSimplex & simplexByIndex(int ind) {return m_simplices[m_indexToSimplex[ind]];}


    //====================================================
    //====================================================
    //==            Complex Construction Functions      ==
    //====================================================
    //====================================================
    /*! Initializer doesn't need to do anything because there's no simplicial structure at this level */
    void init() {}

    /*! Finalizer needs to build m_simplices from simplexSet and generate the index->simplex map*/
    void finalize();

    /*! Since this is the bottom level all add needs to do is check if the simplex exists and index it*/
    int add(NSimplex & simplex);
    //====================================================
    //====================================================
    //==            Geometric Functions                 ==
    //====================================================
    //====================================================
    /*! 0 forms all have volume 1*/
    void computeVolume(NSimplex & s) {s.volume = 1;}



    /*! Compute the volume of the dual polytope by adding the volume of the disjoint components of the dual polytope that are passed down from above*/
    void genDualVolume(NSimplex & s, std::vector<Vector> & clist);
    void genDualVolume(NSimplex & s, std::array<Vector,DimTraits::Top+1> & clist);

    /*! The center of a 0 form is its vertex*/
    void computeCircumcenter(NSimplex & s) {s.center = m_vertices[s[0]];}


    /*! This function should only be called by a higher order simplicial complex, one which will be simply passing down their boundary state which should be blindly followed*/
    void setInterior(int index) {m_interior.diagonal()(index)=0;}

    std::vector<NSimplex > m_simplices;
    DiagonalMatrix m_interior;
    std::vector<Vector> m_vertices;
    std::set<NSimplex> m_simplexSet;
    std::map<int,int> m_indexToSimplex;

};




template <typename NT, typename DT> auto SimplicialComplexPrivateBase<NT,DT>
::add(NSimplex & simplex) -> int {
    simplex.setIndex(-1);
    typename std::set<NSimplex>::const_iterator it = m_simplexSet.find(simplex);
    if(it == m_simplexSet.end()) {
        simplex.setIndex(m_simplexSet.size());
        //m_simplices.push_back(simplex);
        m_simplexSet.insert(simplex);
        return simplex.Index();
    } else {
        return it->Index();
    }
}

template <typename NT, typename DT> auto SimplicialComplexPrivateBase<NT,DT>
::finalize() -> void {
    m_simplices.resize(m_simplexSet.size());
    m_interior.resize(m_simplexSet.size());
    m_interior.setIdentity();
    std::copy(m_simplexSet.begin(), m_simplexSet.end(), m_simplices.begin());
    for(auto&& s: m_simplices)
    {
        computeVolume(s);
        computeCircumcenter(s);
    }
    for(int i=0; i < m_simplices.size(); ++i) {
        m_indexToSimplex[m_simplices[i].Index()] = i;
    }
}





template <typename NT, typename DT> auto SimplicialComplexPrivateBase<NT,DT>
::genDualVolume(NSimplex & s, std::array<Vector,DimTraits::Top+1> & clist) -> void {
    clist[0] = s.center;
    Eigen::Matrix<Scalar,EmbeddedDim,DimTraits::Top> m;
    auto&& origin = clist.back();
    for(int i=0; i < DimTraits::Top; ++i)
    {
        m.col(i) = this->m_vertices[i] - origin;
    }
    s.dualVolume += mtao::volume(m);
}
































/*!  Recursive case for SimplicialComplex class
  * Builds the data necessary for just this level of simplicial complex and then tosses data to the lower complex.
  * There's some necessary intermingling with lower complexes to build certain quantities  like dual volumes which aren't purely top-down quantities
  */
template <typename NT, typename DT>
class SimplicialComplexPrivate: public
        std::conditional<!std::is_same<DT,typename DT::BottomTraits::UpperTraits>::value
        ,SimplicialComplexPrivate<NT,typename DT::LowerTraits>
        , SimplicialComplexPrivateBase<NT,typename DT::BottomTraits>
        >::type
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW
    template <int M>
    struct TraitsContainer{
        typedef mtao_internal::dimensional_traits<DT::Top,M> dimtype;
        typedef typename std::conditional<M==0
                , SimplicialComplexPrivateBase<NT,typename DT::BottomTraits >
                , SimplicialComplexPrivate<NT,dimtype >
                >::type complextype;
        typedef Simplex<NT,dimtype> simplextype;
    };
    typedef NT NumTraits;
    typedef DT DimTraits;
    static const int Dim = DimTraits::Dim;
    static const int N = DimTraits::Dim;
    static const int TopDim = DimTraits::Top;
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    typedef typename TraitsContainer<Dim>::simplextype NSimplex;
    typedef typename TraitsContainer<N-1>::complextype SCm1;
    typedef typename TraitsContainer<0>::complextype SC0;
    typedef Eigen::Matrix<Scalar,EmbeddedDim, N+1> WhitneyBasis;
    typedef Eigen::Matrix<Scalar,N+1, 1> BarycentricCoordinates;
    typedef Eigen::Matrix<Scalar,N+1, 1> WhitneyCoefficients;
protected:
    //====================================================
    //====================================================
    //==               Constructors                     ==
    //====================================================
    //====================================================
    /*! Empty Constructor*/
    SimplicialComplexPrivate() {}
    SimplicialComplexPrivate(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices);
    SimplicialComplexPrivate(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices);
    SimplicialComplexPrivate(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices);
    SimplicialComplexPrivate(const std::vector<std::array<unsigned int, Dim+1> > & tuples, const std::vector<Vector> & vertices);




    //====================================================
    //====================================================
    //==            Complex Construction Functions      ==
    //====================================================
    //====================================================
    void init();
    void finalize();
    int add(NSimplex & simplex);
    NSimplex & simplexByIndex(int ind) {return m_simplices[m_indexToSimplex.at(ind)];}
    const NSimplex & simplexByIndex(int ind) const {return m_simplices[m_indexToSimplex.at(ind)];}
    std::map<int,int> m_indexToSimplex;




    //====================================================
    //====================================================
    //==            Geometric Functions                 ==
    //====================================================
    //====================================================
    Eigen::Matrix<Scalar,EmbeddedDim,N> simplexToBaryMatrix(const NSimplex & s) const;
    void computeCircumcenter(NSimplex & s);
    //template <bool Signed = (N == EmbeddedDim)> void computeVolume(NSimplex & s);
    void genDualVolume(NSimplex & s, std::array<Vector,DimTraits::Top+1> & clist);

    template <bool Signed = (N == EmbeddedDim)>
    void computeVolume(NSimplex & s) {
        auto&& m = simplexToBaryMatrix(s);
        s.volume = mtao::volume(m);
    }

    //====================================================
    //====================================================
    //==            Topological Functions               ==
    //====================================================
    //====================================================
    void createBoundary(NSimplex & simplex);
    void setInterior(int index);

    //====================================================
    //====================================================
    //==            Local Accessors                     ==
    //====================================================
    //====================================================

    Eigen::Matrix<Scalar, EmbeddedDim, Dim+1> verticesBySimplex(const NSimplex & s) const;

    //====================================================
    //====================================================
    //==            Interpolation                       ==
    //====================================================
    //====================================================

    void genWhitneyBases();
    BarycentricCoordinates barycentricCoords(const NSimplex & s, const Vector & v) const;
    bool inSimplex(const NSimplex & s, const Vector & v) const;
    void projectToSimplexInPlace(const NSimplex & s, Vector & v) const;

    std::vector<NSimplex > m_simplices;
    std::vector<Triplet > m_boundaryTriplets;
    //    std::vector< std::array<Vector,N+1 > > m_whitneyBases;
    std::vector< WhitneyBasis > m_whitneyBases;
    std::vector< WhitneyBasis > m_whitneyCenters;
    SparseMatrixColMajor m_boundary;//Rows are n-1 simplices cols are n simplices
    DiagonalMatrix m_interior;
    std::set<NSimplex> m_simplexSet;






};









//====================================================
//====================================================
//==               Constructors                     ==
//====================================================
//====================================================

/*! Constructor that takes in a vector of indexsets and a vector of vertices, should be the standard way to build simplicial complexes
 * \param tuples A collection of n-tuples that form each simplex
 * \param vertices The set of vertices that form the simplicial complex
 */
template <typename NT, typename DT> SimplicialComplexPrivate<NT,DT>::
SimplicialComplexPrivate(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices) {
    m_simplices = simplices;
    SC0::m_vertices=vertices;
    init();
}

template <typename NT, typename DT> SimplicialComplexPrivate<NT,DT>::
SimplicialComplexPrivate(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices) {
    SC0::m_vertices=vertices;
    m_simplices.resize(tuples.size());
    std::transform(tuples.begin(), tuples.end(), m_simplices.begin(),
                   [](const mtao::IndexSet<N+1> & is) -> NSimplex
    {
        return NSimplex(is);
    });

    init();
}

template <typename NT, typename DT> SimplicialComplexPrivate<NT,DT>::
SimplicialComplexPrivate(const std::vector<std::array<unsigned int, Dim+1> > & tuples, const std::vector<Vector> & vertices) {
    SC0::m_vertices=vertices;
    m_simplices.resize(tuples.size());
    std::transform(tuples.begin(), tuples.end(), m_simplices.begin(),
                   [](const std::array<unsigned int, N+1> & is) -> NSimplex
    {
        return NSimplex(static_cast<const mtao::IndexSet<N+1> >(is));
    });

    init();
}

/*! Constructor that takes an unstructured vector of indices and vertices, assuming every N+1 indices form a simplex
 * \param tuples A vector of indices ordered such that [(N+1)*i:(N+1)*(i+1)] is a simplex for each i
 * \param vertices The set of vertices that form the simplicial complex
 */
template <typename NT, typename DT> SimplicialComplexPrivate<NT,DT>::
SimplicialComplexPrivate(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices) {
    SC0::m_vertices=vertices;
    //assert(tuples.size() % N+1 == 0, "Input tuples are incorrectly sized (indices % dimension+1 != 0)");
    m_simplices.resize(tuples.size()/(N+1));
    for(int i=0; i < m_simplices.size(); ++i)
    {
        for(int j=0; j < N+1; ++j){
            m_simplices[i][j] = tuples[(N+1)*i + j];
        }
        //        std::copy(&tuples[(N+1)*i], &tuples[(i+1)*(N+1)], m_simplices[i].getIndexSet().begin());
    }
    m_simplices.resize(tuples.size());
    std::copy(tuples.begin(), tuples.end(), m_simplices.begin());

    init();
}



//====================================================
//====================================================
//==            Complex Construction Functions      ==
//====================================================
//====================================================


template <typename NT, typename DT>
void SimplicialComplexPrivate<NT,DT>::init()
{
    int index=0;
    std::sort(m_simplices.begin(), m_simplices.end());
    //std::unique leaves the iterator at the end to do erase to the end
    m_simplices.erase(std::unique(m_simplices.begin(), m_simplices.end()), m_simplices.end());

    decltype(m_simplices)(m_simplices).swap(m_simplices);


    for(auto && simplex: m_simplices)
    {
        simplex.setIndex(index++);
        createBoundary(simplex);
    }

    finalize();
    //This depends on having the boundary structure, but should be run
    //top down and only once due to the maintenance of circumcenter list
    //so it doesn't fit in head-recursive finalize
    //std::vector<Vector> clist(N+1,Vector::Zero());
    std::array<Vector,DimTraits::Top+1> clist;
    for(auto&& s: m_simplices)
    {
        clist.back() = s.Center();
        s.dualVolume = 1;
        for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it) {
            SCm1::genDualVolume(SCm1::simplexByIndex(it.row()), clist);
        }
        //genDualVolume(s, clist);
    }
    //Extraneous computation to let the DEC side zero things out easier
    SparseMatrixColMajor B = m_boundary.transpose();//make it so the innervectors are the top simplices
    m_interior.resize(m_simplices.size());
    m_interior.setIdentity();
    for(int i = 0; i < B.outerSize(); ++i) {
        if(B.innerVector(i).nonZeros() == 1) {
            SCm1::setInterior(i);
            /*
            m_interior.diagonal()(
                    typename decltype(B)::InnerIterator(B,i).row()
                    )=0;
                    */
        }
    }

}

template <typename NT, typename DT>
void SimplicialComplexPrivate<NT,DT>::finalize()
{

    SCm1::finalize();//need to finalize before boundary, which depends on the smaller pieces

    if(N < TopDim) {
        m_simplices.resize(m_simplexSet.size());
        std::copy(m_simplexSet.begin(), m_simplexSet.end(), m_simplices.begin());
    }
    m_simplexSet.clear();
    std::sort(m_simplices.begin(), m_simplices.end());

    for(auto&& s: m_simplices) {
        computeVolume(s);
        computeCircumcenter(s);
    }


    for(int i=0; i < m_simplices.size(); ++i) {
        m_indexToSimplex[m_simplices[i].Index()] = i;
    }

    m_boundary.resize(SCm1::m_simplices.size(), m_simplices.size());
    m_boundary.setFromTriplets(m_boundaryTriplets.begin(), m_boundaryTriplets.end());
    m_interior.resize(m_simplices.size());
    m_interior.setIdentity();
    genWhitneyBases();



}

//Check whether this simplex already exists
//if it does exist set the index so the owner of the simplex knows where it belongs
//if it doesn't exist push it into the list
template <typename NT,typename DT> auto SimplicialComplexPrivate<NT,DT>
::add(NSimplex & simplex) -> int
{
    simplex.setIndex(-1);
    typename std::set<NSimplex>::const_iterator it = m_simplexSet.find(simplex);
    if(it == m_simplexSet.end())
    {
        simplex.setIndex(m_simplexSet.size());
        //m_simplices.push_back(simplex);
        createBoundary(simplex);
        m_simplexSet.insert(simplex);
        return simplex.Index();
    }
    else
    {
        return it->Index();
    }


}




//====================================================
//====================================================
//==            Geometric Functions                 ==
//====================================================
//====================================================

template <typename NT, typename DT> auto SimplicialComplexPrivate<NT,DT>
::simplexToBaryMatrix(const NSimplex & s) const -> Eigen::Matrix<Scalar,EmbeddedDim,N>
{
    Eigen::Matrix<Scalar,EmbeddedDim,N> m;
    auto&& origin = this->m_vertices[s[N]];
    for(int i=0; i < N; ++i) {
        m.col(i) = this->m_vertices[s[i]] - origin;
    }
    return m;
}

template <typename NT, typename DT> auto SimplicialComplexPrivate<NT,DT>
::computeCircumcenter(NSimplex & s) -> void
{
    //2<V,C> + r = sum(V.array()^2,by col).transpose()
    //switch to homo coords => bottomw ones and right ones, zero bottom left
    auto&& m = simplexToBaryMatrix(s);
    Eigen::Matrix<Scalar,N+1,N+1> A = Eigen::Matrix<Scalar,N+1,N+1>::Ones();

    A.topLeftCorner(N,N) = 2*(m.transpose()*m).eval();
    A.coeffRef(N,N) = 0;
    Eigen::Matrix<Scalar,N+1,1> b = Eigen::Matrix<Scalar,N+1,1>::Ones();
    b.topRows(N) = m.array().square().colwise().sum().transpose();
    A.ldlt().solveInPlace(b);
    s.center = this->m_vertices[s[N]] + m*b.topRows(N);

}

template <typename NT, typename DT> auto SimplicialComplexPrivate<NT,DT>
::genDualVolume(NSimplex & s, std::array<Vector,DimTraits::Top+1> & clist) -> void
{
    clist[N] = s.center;
    Eigen::Matrix<Scalar,EmbeddedDim,DimTraits::Top-N> m;
    auto&& origin = clist.back();
    for(int i=N; i < DimTraits::Top; ++i) {
        m.col(i-N) = clist[i] - origin;
        s.dualVolume += mtao::volume(m);
    }
    for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it) {
        SCm1::genDualVolume(SCm1::simplexByIndex(it.row()), clist);
    }
}


//====================================================
//====================================================
//==            Topological Functions               ==
//====================================================
//====================================================

template <typename NT, typename DT> auto SimplicialComplexPrivate<NT,DT>
::setInterior(int index) -> void {
    m_interior.diagonal()(index)=0;
    for(typename decltype(m_boundary)::InnerIterator it(m_boundary,index);it;++it)
    {
        //if it hasn't been set already
        if(SCm1::m_interior.diagonal()(it.row()) == 1)
            SCm1::setInterior(it.row());
    }
}



template <typename NT,typename DT> auto SimplicialComplexPrivate<NT,DT>
::createBoundary(NSimplex& simplex) -> void
{

    mtao::IndexSet<N> index;

    for ( int i = 0; i <= N ;  ++i) {
        //set a new set of indices
        int k=0;
        for(int j=0; j<N; ++j,++k) {
            index[k]=simplex[(i+j)%(N+1)];
        }
        Simplex<NumTraits,typename DT::LowerTraits > target(index);
        //add the simplex created by add

        if(N != 1) {
            m_boundaryTriplets.push_back(Triplet(
                                             SCm1::add(target),
                                             simplex.Index(),
                                             simplex.isSameSign(target)?1:-1));
        }
        else
        {
            m_boundaryTriplets.push_back(Triplet(
                                             SCm1::add(target),
                                             simplex.Index(),
                                             (i==1)?1:-1));
        }

    }
}


//====================================================
//====================================================
//==            Local Accessors                     ==
//====================================================
//====================================================

template <typename NT,typename DT> auto SimplicialComplexPrivate<NT,DT>
::verticesBySimplex(const NSimplex & s) const
-> Eigen::Matrix<Scalar, EmbeddedDim, Dim+1> {
    Eigen::Matrix<Scalar, EmbeddedDim, Dim+1> m;
    for(int i=0; i <= Dim;++i) {
        m.col(i) = SC0::m_vertices[s[i]];
    }
    return m;
}


//====================================================
//====================================================
//==            Interpolation                       ==
//====================================================
//====================================================

template <typename NT, typename DT>
void SimplicialComplexPrivate<NT,DT>::genWhitneyBases()
{
    m_whitneyBases.resize(m_simplices.size());
    m_whitneyCenters.resize(m_simplices.size());
    //std::cout << "DIM: " << Dim << std::endl;
    for(auto&& s: m_simplices) {
        auto&& basis = m_whitneyBases[s.Index()];
        auto&& centers = m_whitneyCenters[s.Index()];
        int sind=0;
        //std::cout << "======" << s.Volume() << std::endl;
        for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it, ++sind) {
            //Assume that the center is circumcenter, which helps build whitney forms
            auto&& sm1 = SCm1::simplexByIndex(it.row());
            //compute the dimension that is missing

            centers.col(sind) = sm1.Center();
            //std::cout << s.Index() << ": " << sind << ") " << centers.col(sind).transpose() << std::endl;//sm1.Center().transpose() << std::endl;
            //basis vector measured from center to end
            basis.col(sind) = (
                        s.Center()-sm1.Center()).normalized();
            basis.col(sind) *= sm1.Volume() / (N * s.Volume());//need to doublecheck that its n-1
        }
        //std::cout << "vvvvvvv\n";
        //std::cout << m_whitneyBases[s.Index()] << std::endl;
    }
}
template <typename NT, typename DT>
auto SimplicialComplexPrivate<NT,DT>::barycentricCoords(const NSimplex & s, const Vector & v) const -> BarycentricCoordinates
{
    //BarycentricCoordinates coeffs;
    auto&& basis = m_whitneyBases[s.Index()];
    auto&& centers = m_whitneyCenters[s.Index()];
    WhitneyBasis m = v.rowwise().replicate(centers.cols()) - centers;
    return m.cwiseProduct(basis).colwise().sum().transpose();




    /*
    int sind=0;
        WhitneyBasis m = verticesBySimplex(s) - v.rowwise().replicate(m.cols());

    for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it, ++sind) {

        coeffs(sind) = m_whitneyBases[s.Index()].col(sind).dot(v - SCm1::simplexByIndex(it.row()).Center());
    }
    */
    //return coeffs;
}

template <typename NT, typename DT>
auto SimplicialComplexPrivate<NT,DT>::inSimplex(const NSimplex & s, const Vector & v) const -> bool
{
    return !( this->barycentricCoords(s,v).array() < 0 ).any();
}
template <typename NT, typename DT>
void SimplicialComplexPrivate<NT,DT>::projectToSimplexInPlace(const NSimplex & s, Vector & v) const
{
    BarycentricCoordinates bary = barycentricCoords(s,v);
    v = Vector::Zero();
    int i = 0;
    for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it, ++i) {
        auto&& sm1 = SCm1::simplexByIndex(it.row());
        v += SC0::m_vertices[s.oppositeIndex(sm1)] * bary(i);//[i]] * bary(i);
    }
}


};



























template <typename NT, int N_>
class SimplicialComplex: public mtao_internal::SimplicialComplexPrivate<NT,mtao_internal::dimensional_traits<N_,N_> > {
private: typedef mtao_internal::SimplicialComplexPrivate<NT,mtao_internal::dimensional_traits<N_,N_> > PrivateParent;
public:
    typedef NT NumTraits;
    typedef mtao_internal::dimensional_traits<N_,N_> DimTraits;
    typedef mtao_internal::SimplicialComplexTraits<NumTraits, DimTraits> SCTraits;
    static const int Dim = N_;
    static const int N = N_;
    static const int TopDim = N_;
    static const int EmbeddedDim = NT::Dim;
    PLAIN_CLASS_NUM_DEFS
    typedef  Simplex<NumTraits, DimTraits> NSimplex;

    template <int M=N>
    struct TraitsContainer{
        typedef mtao_internal::dimensional_traits<DimTraits::Top,M> dimtype;
        typedef typename std::conditional<M==0
        , mtao_internal::SimplicialComplexPrivateBase<NT,typename DimTraits::BottomTraits >
        , mtao_internal::SimplicialComplexPrivate<NT,dimtype >
        >::type complextype;
        typedef Simplex<NT,dimtype> simplextype;

    };
    //====================================================
    //====================================================
    //==               Constructors                     ==
    //====================================================
    //====================================================
    SimplicialComplex() {}
    SimplicialComplex(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices)
        : PrivateParent(simplices,vertices)
    {
    }
    SimplicialComplex(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }
    SimplicialComplex(const std::vector<std::array<unsigned int, N+1> > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }


    SimplicialComplex(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices)
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


typedef SimplicialComplex<mtao_internal::num_traits<float, 3>, 2> TriangleMeshf;
typedef SimplicialComplex<mtao_internal::num_traits<float, 3>, 3> TetrahedralMeshf;
typedef SimplicialComplex<mtao_internal::num_traits<double, 3>, 2> TriangleMeshd;
typedef SimplicialComplex<mtao_internal::num_traits<double, 3>, 3> TetrahedralMeshd;
//Default to double :)
typedef SimplicialComplex<mtao_internal::num_traits<double, 3>, 2> TriangleMesh;
typedef SimplicialComplex<mtao_internal::num_traits<double, 3>, 3> TetrahedralMesh;






#endif
