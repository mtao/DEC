#ifndef _SIMPLICIAL_COMPLEX_H_
#define _SIMPLICIAL_COMPLEX_H_

#include "simplex.hpp"
#include <vector>
#include <set>

//Input N-simplices to generate the whole simplicial complex
//Assumes that the input mesh of simplices is manifold

//NumTraits::N identifies the dimension of the space
//N determines the maximal simplex dimension

template <typename NT, int N>
class SimplicialComplex: public SimplicialComplex<NT,N-1>
{
public:
    typedef NT NumTraits;
    static const unsigned int Dim = N;
    virtual unsigned int topDim() const{return Dim;}
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
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
        std::copy(tuples.begin(), tuples.end(), m_simplices.begin());

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
    //Builds the n-1 simplices and n <-> n-1 simplices relationships
    void init();
    int add(NSimplex & simplex);
    int createBoundary(NSimplex & simplex);
    //const std::vector<NSimplex> & Simplices() const {return m_simplices;}

protected:
    std::vector<NSimplex > m_simplices;
    void finalize();
private:
    typedef SimplicialComplex<NT,N-1> SCm1;
    std::vector<Triplet > m_boundaryTriplets;
    SparseMatrixColMajor m_boundary;//Rows are n-1 simplices cols are n simplices
    std::set<NSimplex> m_simplexSet;

};

//A 0-simplicial complex only manages
template <typename NT>
class SimplicialComplex<NT,0>
{
public:
    typedef NT NumTraits;
    static const unsigned int Dim = 0;
    virtual unsigned int topDim() const {return Dim;}
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::Triplet Triplet;
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

protected:
    void init() {}
    void finalize()
    {
        m_simplices.resize(m_simplexSet.size());
        std::copy(m_simplexSet.begin(), m_simplexSet.end(), m_simplices.begin());
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
protected:
    std::vector<NSimplex > m_simplices;
    std::vector<Vector> m_vertices;
    std::set<NSimplex> m_simplexSet;

};



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

        m_boundaryTriplets.push_back(Triplet(
                                         SimplicialComplex<NT,N-1>::add(target),
                                         simplex.Index(),
                                         simplex.isSameSign(target)?1:-1));
        //std::cout << N << ": "
        //          << Simplices<N-1>()[m_boundaryTriplets.back().row()] << " "
        //          << Simplices<N>()[m_boundaryTriplets.back().col()] << " " << std::endl;


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

}

template <typename NT, int N>
void SimplicialComplex<NT,N>::finalize()
{
    SCm1::finalize();
    if(N < topDim())
    {
    m_simplices.resize(m_simplexSet.size());
    std::copy(m_simplexSet.begin(), m_simplexSet.end(), m_simplices.begin());
    }
    m_simplexSet.clear();
    m_boundary.resize(SCm1::m_simplices.size(), m_simplices.size());
    m_boundary.setFromTriplets(m_boundaryTriplets.begin(), m_boundaryTriplets.end());
    //std::cout << m_boundary << std::endl;
}



typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMesh;








#endif
