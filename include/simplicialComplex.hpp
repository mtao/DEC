#ifndef _SIMPLICIAL_COMPLEX_H_
#define _SIMPLICIAL_COMPLEX_H_

#include "simplex.hpp"
#include <vector>


//Input N-simplices to generate the whole simplicial complex
//Assumes that the input mesh of simplices is manifold

//NumTraits::N identifies the dimension of the space
//N determines the maximal simplex dimension

template <typename NT, int N>
class SimplicialComplex: public SimplicialComplex<NT,N-1>
{
public:
    typedef NT NumTraits;
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
                                      *(reinterpret_cast<mtao::IndexSet<N+1> *>(&tuples[i]))
                                      )
                                  );

        }
        m_simplices.resize(tuples.size());
        std::copy(tuples.begin(), tuples.end(), m_simplices.begin());

        init();
    }
    template <int M>
    const std::vector<Simplex<NumTraits, M> > & Simplices() const
    {
        static_assert( M <= N, "Tried to access a set of simplices of a dimension too high");
        return SimplicialComplex<NT,M>::m_simplices;
    }
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

};



//A 0-simplicial complex only manages
template <typename NT>
class SimplicialComplex<NT,0>
{
public:
    typedef NT NumTraits;
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
    const std::vector<Vector> & Vertices(){return m_vertices;}

protected:
    void init() {}
    void finalize()
    {}
    int add(NSimplex & simplex)
    {
        for(auto it=m_simplices.begin();
            it!= m_simplices.end();++it)

        {
            if(*it==simplex)
            {
                simplex.setIndex(it->Index());
                break;
            }
        }
        if(simplex.Index()==-1)
        {
            simplex.setIndex(m_simplices.size());
            m_simplices.push_back(simplex);

        }
        return simplex.Index();
    }
protected:
    std::vector<NSimplex > m_simplices;
    std::vector<Vector> m_vertices;

};

//Check whether this simplex already exists
//if it does exist set the index so the owner of the simplex knows where it belongs
//if it doesn't exist push it into the list 
template <typename NT,int N>
int SimplicialComplex<NT,N>::add(NSimplex & simplex)
{
    simplex.setIndex(-1);
    for(auto it=m_simplices.begin();
        it!= m_simplices.end();++it)

    {
        if(*it==simplex)
        {
            simplex.setIndex(it->Index());
            break;
        }
    }
    if(simplex.Index()==-1)
    {
        simplex.setIndex(m_simplices.size());
        m_simplices.push_back(simplex);
        createBoundary(simplex);

    }
    return simplex.Index();


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
    std::vector<Triplet > tripList;
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
    m_boundary.resize(SCm1::m_simplices.size(), m_simplices.size());
    m_boundary.setFromTriplets(m_boundaryTriplets.begin(), m_boundaryTriplets.end());
    //std::cout << m_boundary;
    SCm1::finalize();
}



typedef SimplicialComplex<NumericalTraits<double, 3>, 2> TriangleMesh;








#endif
