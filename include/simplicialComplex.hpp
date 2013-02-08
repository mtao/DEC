#ifndef _SIMPLICIAL_COMPLEX_H_
#define _SIMPLICIAL_COMPLEX_H_

#include "simplex.hpp"
#include <vector>
#include <array>
#include <algorithm>
#include <set>

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
template <typename NT, int N>
class SimplicialComplex: public SimplicialComplex<NT,N-1>
{
public:
    typedef NT NumTraits;
    static const int Dim = N;
    virtual int TopDim() const { return Dim; }
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef typename NumTraits::SparseMatrixColMajor SparseMatrixColMajor;
    typedef  Simplex<NumTraits, N> NSimplex;
    SimplicialComplex() {}
    SimplicialComplex(const std::vector<NSimplex > & simplices, const std::vector<Vector> & vertices)
    {
        m_simplices = simplices;
        SimplicialComplex<NT,0>::m_vertices=vertices;
        init();
    }
    //N+1 vertices on an N-simplex
    SimplicialComplex(const std::vector<mtao::IndexSet<N+1> > & tuples, const std::vector<Vector> & vertices)
    {
        SimplicialComplex<NT,0>::m_vertices=vertices;
        m_simplices.resize(tuples.size());
        std::transform(tuples.begin(), tuples.end(), m_simplices.begin(),
                [](const mtao::IndexSet<N+1> & is) -> NSimplex
                {
                return NSimplex(is);
                });

        init();
    }

    SimplicialComplex(const std::vector<unsigned int > & tuples, const std::vector<Vector> & vertices)
    {
        SimplicialComplex<NT,0>::m_vertices=vertices;
        //assert(tuples.size() % N+1 == 0, "Input tuples are incorrectly sized (indices % dimension+1 != 0)");
        m_simplices.resize(0);
        m_simplices.reserve(tuples.size()/(N+1));
        for(int i=0; i < tuples.size(); i+=N+1)
        {
            m_simplices.push_back(Simplex<NumTraits, N>
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
    std::vector<unsigned int> simplicesToArray() const

    {
        static_assert(M > 0 && M <= Dim, "Inappropriate dims asserted");
        auto&& simplices = this->template constSimplices<M>();
        std::vector<unsigned int> ret(simplices.size() * (M+1));
        for(auto&& s: simplices)
        {
            std::copy(s.getIndexSet().cbegin(),s.getIndexSet().cend(), ret.begin()+s.Index() * (M+1));
            if(s.isNegative())
            {
                unsigned int tmp =  ret[s.Index() * (M+1)+1];
                ret[s.Index() * (M+1)+1] =  ret[s.Index() * (M+1)];
                ret[s.Index() * (M+1)] =  tmp;
            }
        }
        return ret;
    }

    template <int M=Dim>
    const SparseMatrix & b() const{return SimplicialComplex<NT,M>::m_boundary;}
    template <int M=Dim>
    const std::vector<Simplex<NumTraits, M> > & constSimplices() const
    {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return SimplicialComplex<NT,M>::m_simplices;
    }
    template <int M=Dim>
    std::vector<Simplex<NumTraits, M> > & simplices()
    {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return SimplicialComplex<NT,M>::m_simplices;
    }
    template <int M=Dim>
    size_t numSimplices()const{return SimplicialComplex<NT,M>::m_simplices.size();}

    template <int M=Dim>
    const DiagonalMatrix & interior() const 
    {return SimplicialComplex<NT,M>::m_interior;}
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
                m.col(i-N) = this->m_vertices[i] - origin;
            }
            s.dualVolume += std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(N);
        }
        for(typename decltype(m_boundary)::InnerIterator it(m_boundary, s.Index()); it; ++it)
        {
            SCm1::genDualVolume(simplices<N-1>()[it.row()], clist);
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
    typedef SimplicialComplex<NT,N-1> SCm1;
    std::vector<Triplet > m_boundaryTriplets;
    std::vector< std::array<Vector,N+1 > > m_whitneyBases;
    SparseMatrixColMajor m_boundary;//Rows are n-1 simplices cols are n simplices
    DiagonalMatrix m_interior;
    std::set<NSimplex> m_simplexSet;



};

//A 0-simplicial complex only manages
//the vertices and some trivalish set of simplices.  vertices without any higher order
//information are considered bad data and though stored here, they do not have simplices 
//associated with them.
template <typename NT>
class SimplicialComplex<NT,0>
{
public:
    typedef NT NumTraits;
    static const int Dim = 0;
    static const int EmbeddedDim = NT::Dim;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::DiagonalMatrix DiagonalMatrix;
    typedef Simplex<NumTraits,0> NSimplex;


    SimplicialComplex() {}
    SimplicialComplex(const std::vector<Vector> & vertices)
    {
        m_vertices=vertices;
    }

    SimplicialComplex(const std::vector<NSimplex >& simplices, const std::vector<Vector> & vertices)
    {
        m_simplices = simplices;
        m_vertices=vertices;
    }
    std::vector<Vector> & vertices(){return m_vertices;}
    const std::vector<Vector> & constVertices()const{return m_vertices;}
    template <int M=0>
    const DiagonalMatrix & interior() const 
    {return SimplicialComplex<NT,M>::m_interior;}

protected:
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
        s.center = m_vertices[s.Index()];
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
        s.dualVolume += std::sqrt((m.transpose()*m).determinant())/mtao::cefactorial(0);
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

};



template <typename NT, int N>
    void SimplicialComplex<NT,N>::genWhitneyBases()
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
                auto&& v = this->vertices()[n];
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
template <typename NT, int N>
auto SimplicialComplex<NT,N>::barycentricCoords(const NSimplex & s, const Vector & v) -> std::array<Scalar,N+1>
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
    template <typename NT,int N>
int SimplicialComplex<NT,N>::add(NSimplex & simplex)
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
    template <typename NT,int N>
int SimplicialComplex<NT,N>::createBoundary(NSimplex& simplex)
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
        Simplex<NumTraits,N-1> target(index);
        //add the simplex created by add

        if(N != 1) {
            m_boundaryTriplets.push_back(Triplet(
                        SimplicialComplex<NT,N-1>::add(target),
                        simplex.Index(),
                        simplex.isSameSign(target)?1:-1));
        }
        else
        {
            m_boundaryTriplets.push_back(Triplet(
                        SimplicialComplex<NT,N-1>::add(target),
                        simplex.Index(),
                        //(simplex.isSameSign(target)==(i==1))?1:-1));
                (i==1)?1:-1));
        }

    }
    return 0;
}
    template <typename NT, int N>
void SimplicialComplex<NT,N>::init()
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
            m_interior.diagonal()(
                    typename decltype(B)::InnerIterator(B,i).row()
                    )=0;


        }
    }

}

template <typename NT, int N>
void SimplicialComplex<NT,N>::finalize()
{

    SCm1::finalize();//need to finalize before boundary, which depends on the smaller pieces

    if(N < TopDim())
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



    m_boundary.resize(SCm1::m_simplices.size(), m_simplices.size());
    m_boundary.setFromTriplets(m_boundaryTriplets.begin(), m_boundaryTriplets.end());
    m_interior.resize(m_simplices.size());
    m_interior.setIdentity();
    genWhitneyBases();
    


}


typedef SimplicialComplex<NumericalTraits<float, 3>, 2> TriangleMeshf;
typedef SimplicialComplex<NumericalTraits<float, 3>, 3> TetrahedralMeshf;
typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMeshd;
typedef SimplicialComplex<NumericalTraits<double, 3>, 3> TetrahedralMeshd;
//Default to double :)
typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMesh;
typedef SimplicialComplex<NumericalTraits<double, 3>, 3> TetrahedralMesh;








#endif
