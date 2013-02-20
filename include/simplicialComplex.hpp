#ifndef _SIMPLICIAL_COMPLEX_H_
#define _SIMPLICIAL_COMPLEX_H_

#include "simplex.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <set>
#include <iostream>

//Input N-simplices to generate the whole simplicial complex
//Assumes that the input mesh of simplices is manifold

//NumTraits::N identifies the dimension of the space
//N determines the maximal simplex dimension

namespace mtao
{
constexpr int cefactorial(int n)
{
    return n > 0 ? n * cefactorial(n-1):1;
}
};
template <int Top_, int Dim_>
struct DimTraits {
    const static int Top = Top_;
    const static int Dim = Dim_;
    typedef DimTraits<Top,Dim-1> LowerTraits;
    typedef DimTraits<Top,Dim+1> UpperTraits;
    typedef DimTraits<Top,Top> TopTraits;
    typedef DimTraits<Top,0> BottomTraits;
};
//A 0-simplicial complex only manages
//the vertices and some trivalish set of simplices.  vertices without any higher order
//information are considered bad data and though stored here, they do not have simplices
//associated with them.
template <typename NT, typename DT>
class SimplicialComplexPrivateBase
{
public:
    template <int M>
    struct SCParent{
        typedef typename std::enable_if<M==0
        , SimplicialComplexPrivateBase<NT,DimTraits<DT::Top, M> > >::type type;
    };
    typedef NT NumTraits;
    static const int Dim = 0;
    static const int TopDim = DT::Top;
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef Simplex<NumTraits,DT> NSimplex;


    SimplicialComplexPrivateBase() {}
    SimplicialComplexPrivateBase(const std::vector<Vector> & vertices)
    {
        m_vertices=vertices;
    }

    SimplicialComplexPrivateBase(const std::vector<NSimplex >& simplices, const std::vector<Vector> & vertices)
    {
        m_simplices = simplices;
        m_vertices=vertices;
    }
    std::vector<Vector> & vertices(){return m_vertices;}
    const std::vector<Vector> & constVertices()const{return m_vertices;}
    template <int M=0>
    const DiagonalMatrix & interior() const
    {return SimplicialComplexPrivateBase<NT,typename SCParent<M>::type >::m_interior;}
    template <int M=0>
    const std::map<int,int> & indexToSimplex() const {return SCParent<M>::type::m_indexToSimplex;}

public://protected:
    void init() {}
    void finalize()
    {
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
    int add(NSimplex & simplex)
    {
        simplex.setIndex(-1);
        typename std::set<NSimplex>::const_iterator it = m_simplexSet.find(simplex);
        if(it == m_simplexSet.end())
        {
            simplex.setIndex(m_simplexSet.size());
            //m_simplices.push_back(simplex);
            m_simplexSet.insert(simplex);
            return simplex.Index();
        }
        else
        {
            return it->Index();
        }
    }
    void computeVolume(NSimplex & s)
    {
        s.volume = 1;
    }
    void computeCircumcenter(NSimplex & s)
    {
        s.center = m_vertices[s[0]];
    }
    template <int M=0>
    typename SCParent<M>::type::NSimplex & simplexByIndex(int ind) {
        return SCParent<M>::type::m_simplices[
            SCParent<M>::type::m_indexToSimplex[ind]
            ];
    }

    void genDualVolume(NSimplex & s, std::vector<Vector> & clist)
    {
        clist[0] = s.center;
        int M = EmbeddedDim;
        typename NumTraits::DynamicMatrix m(M,clist.size()-1);
        auto&& origin = clist.back();
        for(int i=0; i < static_cast<int>(clist.size())-1; ++i)
        {
            m.col(i) = this->m_vertices[i] - origin;
        }
        s.dualVolume += std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(TopDim-1);
    }
    void setInterior(int index)
    {
        m_interior.diagonal()(index)=0;
    }

protected:
    std::vector<NSimplex > m_simplices;
    DiagonalMatrix m_interior;
    std::vector<Vector> m_vertices;
    std::set<NSimplex> m_simplexSet;
    std::map<int,int> m_indexToSimplex;

};
template <typename NT, typename DT>
class SimplicialComplexPrivate: public
        std::conditional<!std::is_same<DT,typename DT::BottomTraits::UpperTraits>::value
        ,SimplicialComplexPrivate<NT,typename DT::LowerTraits>
        , SimplicialComplexPrivateBase<NT,typename DT::BottomTraits>
        >::type
{
public:
    template <int M>
    struct SCParent{
        typedef typename std::conditional<M==0
        , SimplicialComplexPrivateBase<NT,DimTraits<DT::Top, M> >
        , SimplicialComplexPrivate<NT,DimTraits<DT::Top, M> >
        >::type type;
    };
    template <int M>
    struct LocalDT{
        typedef DimTraits<DT::Top,M> type;
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
    typedef  Simplex<NumTraits, DT> NSimplex;
    SimplicialComplexPrivate() {}
    SimplicialComplexPrivate(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices)
    {
        m_simplices = simplices;
        SC0::m_vertices=vertices;
        init();
    }
    //N+1 vertices on an N-simplex
    SimplicialComplexPrivate(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
    {
        SC0::m_vertices=vertices;
        m_simplices.resize(tuples.size());
        std::transform(tuples.begin(), tuples.end(), m_simplices.begin(),
                [](const mtao::IndexSet<N+1> & is) -> NSimplex
                {
                return NSimplex(is);
                });

        init();
    }

    SimplicialComplexPrivate(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices)
    {
        SC0::m_vertices=vertices;
        //assert(tuples.size() % N+1 == 0, "Input tuples are incorrectly sized (indices % dimension+1 != 0)");
        m_simplices.resize(0);
        m_simplices.reserve(tuples.size()/(N+1));
        for(int i=0; i < tuples.size(); i+=N+1)
        {
            m_simplices.push_back(NSimplex
                                  (
                                      reinterpret_cast<mtao::IndexSet<N+1> >(&tuples[i])
                                      )
                                  );

        }
        m_simplices.resize(tuples.size());
        std::copy(tuples.begin(), tuples.end(), m_simplices.begin());

        init();
    }


    template <int M=Dim>
    const SparseMatrix & b() const{return SCParent<M>::type::m_boundary;}
    template <int M=Dim>
    const std::vector<Simplex<NumTraits, typename LocalDT<M>::type > > & constSimplices() const
    {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return SCParent<M>::type::m_simplices;
    }
    template <int M=Dim>
    std::vector<Simplex<NumTraits, typename LocalDT<M>::type > > & simplices()
    {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return SCParent<M>::type::m_simplices;
    }
    template <int M=Dim>
    size_t numSimplices()const{return SCParent<M>::type::m_simplices.size();}

    template <int M=Dim>
    const DiagonalMatrix & interior() const 
    {return SCParent<M>::type::m_interior;}

    template <int M=Dim>
    typename SCParent<M>::type::NSimplex & simplexByIndex(int ind) {
        return SCParent<M>::type::m_simplices[
            SCParent<M>::type::m_indexToSimplex[ind]
            ];
    }
    template <int M=Dim>
    const std::map<int,int> & indexToSimplex() const {return SCParent<M>::type::m_indexToSimplex;}
protected:
    //Builds the n-1 simplices and n <-> n-1 simplices relationships
    void init();
    int add(NSimplex & simplex);
    int createBoundary(NSimplex & simplex);
    void genWhitneyBases();
    std::array<Scalar,N+1> barycentricCoords(const NSimplex & s, const Vector & v);

    Eigen::Matrix<Scalar,EmbeddedDim,N> simplexToBaryMatrix(const mtao::IndexSet<N+1> & s)
    {
        Eigen::Matrix<Scalar,EmbeddedDim,N> m;
        auto&& origin = this->m_vertices[s[N]];
        for(int i=0; i < N; ++i)
        {
            m.col(i) = this->m_vertices[s[i]] - origin;
        }
        return m;
    }
    Eigen::Matrix<Scalar,EmbeddedDim,N> simplexToBaryMatrix(NSimplex & s)
    {
        return simplexToBaryMatrix(s.getIndexSet());
    }
    template <bool Signed = (N == EmbeddedDim)>
    Scalar computeVolume(const mtao::IndexSet<N+1> & s)
    {
        auto&& m = simplexToBaryMatrix(s);
        if(Signed)
        {
            return m.determinant()/mtao::cefactorial(N);
        }
        else
        {
            return std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(N);
        }
    }

    template <bool Signed = (N == EmbeddedDim)>
    void computeVolume(NSimplex & s)
    {
        s.volume = computeVolume(s.getIndexSet());
    }

    void computeCircumcenter(NSimplex & s)
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


    void genDualVolume(NSimplex & s, std::vector<Vector> & clist)
    {
        /*
        for(int i=0; i < TopDim - Dim; ++i)
        std::cout << " ";
        std::cout << s.Index() << " " << s.center.transpose() << std::endl;
        */
        clist[N] = s.center;

        if(N == clist.size()-1)
            s.dualVolume = 1;
        else
        {
            int M = EmbeddedDim;
            typename NumTraits::DynamicMatrix m(M,clist.size()-N-1);
            auto&& origin = clist.back();
            for(int i=N; i < static_cast<int>(clist.size())-1; ++i)
            {
                m.col(i-N) = clist[i] - origin;
                /*
                std::cout << (SC0::m_vertices[i] - origin).eval().transpose() <<  "======";
                std::cout << SC0::m_vertices[i].transpose() << "----" << origin.transpose() <<std::endl;
                */
            }
            /*
            if(N==1 && false)
           {
                for(int i=N; i < clist.size(); ++i) {
        for(int i=0; i < clist.size(); ++i)
            std::cout << " ";
        std::cout << i << "---- " << clist[i].transpose() << std::endl;
                }

        for(int i=0; i < TopDim - Dim; ++i)
        std::cout << " ";
                std::cout <<"Resulting mat: "  << m.transpose() << std::endl;
        for(int i=0; i < TopDim - Dim; ++i)
        std::cout << " ";
            std::cout << "Volume computed:" << s.Index()  << " " << std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(TopDim - N) << std::endl;
            }
            */
            s.dualVolume += std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(TopDim - N);
        }
        for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it)
        {
            SCm1::genDualVolume(simplexByIndex<N-1>(it.row()), clist);
        }
    }


    /*
    void computeDualVolume(NSimplex & s, mtao::IndexSet<EmbeddedDim> & ind)
    {

    }
    */
    void setInterior(int index)
    {
        m_interior.diagonal()(index)=0;
        for(typename decltype(m_boundary)::InnerIterator it(m_boundary,index);it;++it)
        {
            //if it hasn't been set already
            if(SCm1::m_interior.diagonal()(it.row()) == 1)
                SCm1::setInterior(it.row());
        }
    }

protected:
    std::vector<NSimplex > m_simplices;
    void finalize();
protected:
    typedef typename SCParent<N-1>::type SCm1;
    typedef typename SCParent<0>::type SC0;
    std::vector<Triplet > m_boundaryTriplets;
    std::vector< std::array<Vector,N+1 > > m_whitneyBases;
    SparseMatrixColMajor m_boundary;//Rows are n-1 simplices cols are n simplices
    DiagonalMatrix m_interior;
    std::set<NSimplex> m_simplexSet;
    std::map<int,int> m_indexToSimplex;



};




template <typename NT, typename DT>
    void SimplicialComplexPrivate<NT,DT>::genWhitneyBases()
{
    m_whitneyBases.resize(m_simplices.size());
    std::transform(m_simplices.begin(), m_simplices.end(), m_whitneyBases.begin(), [&](const NSimplex & s) 
            //-> std::map<int,Vector>
            -> std::array<Vector,N+1>
        {
            //std::vector<int,Vector> mymap;
            std::array<Vector,N+1> myarr;
            int sind=0;
            for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it, ++sind)
            {
                //Assume that the center is circumcenter, which helps build whitney forms
                auto&& sm1 = simplices<N-1>()[it.row()];
                //compute the dimension that is missing
                int n=s[N];
                for(int i=0; i < N; ++i)
                {
                    if(s[i] != sm1[i])
                    {
                        n = s[i];
                        break;
                    }
                }
                auto&& m1cc = sm1.Center();
                auto&& v = SC0::m_vertices[n];
                Vector d = s.Center()-m1cc;
                Scalar d2n = d.squaredNorm();
                Scalar scale = (v-m1cc).dot(d)/(d2n*d2n);
                /*
                mymap.insert(std::pair<int,Vector>(n,
                        scale*d));
                */
                myarr[sind] = scale * d;


        }
//        return mymap;
        return myarr;
    });
}
template <typename NT, typename DT>
auto SimplicialComplexPrivate<NT,DT>::barycentricCoords(const NSimplex & s, const Vector & v) -> std::array<Scalar,N+1>
{
    std::array<Scalar,N+1> myvec;
    int sind=0;
    for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it, ++sind)
    {
        myvec[sind] = m_whitneyBases[s.Index()][sind].dot(v - simplices<N-1>()[it.row()]);
    }
}
//Check whether this simplex already exists
//if it does exist set the index so the owner of the simplex knows where it belongs
//if it doesn't exist push it into the list
    template <typename NT,typename DT>
int SimplicialComplexPrivate<NT,DT>::add(NSimplex & simplex)
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
    template <typename NT,typename DT>
int SimplicialComplexPrivate<NT,DT>::createBoundary(NSimplex& simplex)
{

    mtao::IndexSet<N> index;

    for ( int i = 0; i <= N ;  ++i)
    {
        //set a new set of indices
        int k=0;
        for(int j=0; j<N; ++j,++k)
        {
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
                        //(simplex.isSameSign(target)==(i==1))?1:-1));
                (i==1)?1:-1));
        }

    }
    return 0;
}
    template <typename NT, typename DT>
void SimplicialComplexPrivate<NT,DT>::init()
{
    int index=0;
    /*
    std::sort(m_simplices.begin(), m_simplices.end(), [](const NSimplex & a, const NSimplex & b) -> bool{
            return a.Index() < b.Index();
            });
            */
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
    std::vector<Vector> clist(N+1,Vector::Zero());
    for(auto&& s: m_simplices)
    {
            genDualVolume(s, clist);
    }
    //Extraneous computation to let the DEC side zero things out easier
    SparseMatrixColMajor B = m_boundary.transpose();//make it so the innervectors are the top simplices
    m_interior.resize(m_simplices.size());
    m_interior.setIdentity();
    for(int i = 0; i < B.outerSize(); ++i)
    {
        if(B.innerVector(i).nonZeros() == 1)
        {
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

    if(N < TopDim)
    {
        m_simplices.resize(m_simplexSet.size());
        std::copy(m_simplexSet.begin(), m_simplexSet.end(), m_simplices.begin());
    }
    m_simplexSet.clear();
    std::sort(m_simplices.begin(), m_simplices.end());

    for(auto&& s: m_simplices)
    {
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
template <typename NT, int N_>
class SimplicialComplex: public SimplicialComplexPrivate<NT,DimTraits<N_,N_> > {
    private: typedef SimplicialComplexPrivate<NT,DimTraits<N_,N_> > PrivateParent;
    public:
    typedef NT NumTraits;
    typedef DimTraits<N_,N_> DimTraits;
    static const int Dim = N_;
    static const int N = N_;
    static const int TopDim = N_;
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    typedef  Simplex<NumTraits, DimTraits> NSimplex;
    SimplicialComplex() {}
    SimplicialComplex(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices)
        : PrivateParent(simplices,vertices)
    {
    }
    SimplicialComplex(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
    }

    SimplicialComplex(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices)
        : PrivateParent(tuples,vertices)
    {
        }

};


typedef SimplicialComplex<NumericalTraits<float, 3>, 2> TriangleMeshf;
typedef SimplicialComplex<NumericalTraits<float, 3>, 3> TetrahedralMeshf;
typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMeshd;
typedef SimplicialComplex<NumericalTraits<double, 3>, 3> TetrahedralMeshd;
//Default to double :)
typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMesh;
typedef SimplicialComplex<NumericalTraits<double, 3>, 3> TetrahedralMesh;








#endif
