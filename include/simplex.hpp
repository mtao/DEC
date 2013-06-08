#ifndef SIMPLEX_H
#define SIMPLEX_H
#include <m_vector>
#include "types.hpp"


namespace DEC{
namespace internal{
template <typename NT, typename DT>
class Complex_<Simplex<NT,DT> >;
template <typename NT, typename DT>
class ComplexBase_<Simplex<NT, DT> >;
};
};

template <typename NT, typename DT>
class Simplex
{
public:
    friend class DEC::internal::Complex_<NT,DT>;
    friend class DEC::internal::ComplexBase_<NT,DT>;
    typedef NT num_traits;
    typedef DT dim_traits;
    enum {dim = DT::dim};
    typedef typename num_traits::Vector Vector;
    typedef typename num_traits::Scalar Scalar;
    bool isNegative( void )  const  { return (m_sign < 0); }
    void negate    ( void )         { m_sign = -m_sign; }
    void resetSign ( void )         { m_sign = 1; }

    inline int Index     ( void )  const  { return m_index; }
    inline void setIndex(int i)              { m_index = i;}
    inline const mtao::IndexSet<Dim+1> & getIndexSet() const {return m_v;}
    Vector center    ( void )  const  { return m_center; }
    Scalar volume    ( void )  const  { return m_volume; }
    Scalar dualVolume( void )  const  { return m_dual_volume; }
    void setCenter    ( void )  const  { return m_center; }
    void setVolume    ( void )  const  { return m_volume; }
    void setDualVolume( void )  const  { return m_dual_volume; }

    /*! Returns the m_index of the m_vertex opposite the simplex*/
    unsigned int oppositeIndex(const Simplex<NT,typename DT::LowerTraits> & below) const {
            for(int i=0; i < Dim; ++i) {
                if(below[i] != m_v[i]) {
                    return m_v[i];
                }
            }
            return m_v[dim];
    }

    Simplex() {}
    Simplex( const mtao::IndexSet<Dim+1> & _v, bool sort = true )
    {
        int i=0;
        m_v = _v;

        if( sort )
        {
            // bubble sort (might as well, N shouldn't be larger than 4)
            int tmp;
            int j;
            for( i = 0; i < Dim; ++i )
            {
                for( j = 0; j < Dim-i; ++j )
                {
                    if( m_v[j+1] < m_v[j] )
                    {
                        tmp = m_v[j];
                        m_v[j] = m_v[j+1];
                        m_v[j+1] = tmp;
                        m_sign *= -1;
                    }
                }
            }
        }

    }

    Simplex( const Simplex& rhs )
    {
        *this = rhs;
    }

    Simplex& operator=( const Simplex& rhs )
    {
        static_assert(false,"Copying simplex is disallowed");
        return (*this);
    }

    unsigned int operator[] (int i)const {
        return m_v[i];
    }
    unsigned int& operator[] (int i) {
        return m_v[i];
    }
    friend inline std::ostream& operator<<(std::ostream& os, const Simplex<NT,DT> & simplex)
    {
        os << "(" << simplex.Index() << ": ";
        if(simplex.isNegative())
        {
            os << "-";
        }
        else
        {
            os << "+";
        }
        os << "{" <<simplex.Volume() << ","<< simplex.DualVolume()<<"}";

        os << simplex.getIndexSet() << ")";
        return os;
    }


protected:
    int m_sign=1;
    size_t m_index=-1;
    Vector center=Vector::Zero();
    Scalar m_volume=0;
    Scalar m_dual_volume=0;
    mtao::IndexSet<Dim+1> m_v;



    /*
     * OPERATORS
     *
     */
public:

    inline bool operator==(  const Simplex<NT,DT>& other ) const
    {
        return this->getIndexSet() == other.getIndexSet();
    }

    inline bool operator!=(  const Simplex<NT,DT>& other ) const
    {
        return (!(*this == other));
    }

    inline bool operator<(  const Simplex<NT,DT>& other ) const
    {
        return this->getIndexSet() < other.getIndexSet();
    }

    inline bool operator>(  const Simplex<NT,DT>& other ) const
    {
        return (other < *this);
    }

    inline bool operator<=(  const Simplex<NT,DT>& other ) const
    {
        return (!(other < this));
    }

    inline bool operator>=(  const Simplex<NT,DT>& other ) const
    {
        return (!(*this < other));
    }

    template <typename DT2>
    inline bool isSameSign( const Simplex<NT,DT2>& other ) const
    {
        return (this->isNegative() == other.isNegative());
    }


};
template<typename NT, typename DT>
inline Simplex<NT,DT> operator-( const Simplex<NT,DT>& simplex )
{
    Simplex<NT,DT> negated( simplex );
    negated.Negate();
    return negated;
}



/*
template <typename NT>
inline Simplex<NT,0> MakeSimplex( int m_v0 )
{
    int m_v[1] = {m_v0};
    return Simplex<NT,0>( m_v );
}

template <typename NT>
inline Simplex<NT,1> MakeSimplex( int m_v0, int m_v1 )
{
    int m_v[2] = {m_v0, m_v1};
    return Simplex<NT,1>( m_v );
}

template <typename NT>
inline Simplex<NT,2> MakeSimplex( int m_v0, int m_v1, int m_v2 )
{
    int m_v[3] = {m_v0, m_v1, m_v2};
    return Simplex<NT,2>( m_v );
}

template <typename NT>
inline Simplex<NT,3> MakeSimplex( int m_v0, int m_v1, int m_v2, int m_v3 )
{
    int m_v[4] = {m_v0, m_v1, m_v2, m_v3};
    return Simplex<NT,3>( m_v, false );
}
*/

#endif
