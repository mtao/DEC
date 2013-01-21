#ifndef SIMPLEX_H
#define SIMPLEX_H
#include <vector>
#include "types.hpp"




template <typename NT, int N>
class SimplicialComplex;

template <typename NT, int DIM=NT::N>
class Simplex
{
public:
    friend class SimplicialComplex<NT,DIM>;
    typedef NT NumTraits;
    static const int Dim = DIM;
    typedef typename NumTraits::Vector Vector;
    typedef typename NumTraits::Scalar Scalar;
    bool IsNegative( void )  const  { return (sign < 0); }
    void Negate    ( void )         { sign = -sign; }
    void ResetSign ( void )         { sign = 1; }

    static int N() {return Dim;}
    inline int Index     ( void )  const  { return index; }
    inline void setIndex(int i)              { index = i;}
    inline const mtao::IndexSet<Dim+1> & getIndexSet() const {return v;}
    Vector Center    ( void )  const  { return center; }
    Scalar Volume    ( void )  const  { return volume; }
    Scalar   DualVolume( void )  const  { return dualVolume; }

    Simplex() {}
    Simplex( const mtao::IndexSet<Dim+1> & _v, bool sort = true )
    {
        int i=0;
        v = _v;

        if( sort )
        {
            // bubble sort (might as well, N shouldn't be larger than 4)
            int tmp;
            int j;
            for( i = 0; i < Dim; ++i )
            {
                for( j = 0; j < Dim-i; ++j )
                {
                    if( v[j+1] < v[j] )
                    {
                        tmp = v[j];
                        v[j] = v[j+1];
                        v[j+1] = tmp;
                        sign *= -1;
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
        sign = rhs.sign;
        v = rhs.v;
        index = rhs.index;
        center = rhs.center;
        volume = rhs.volume;
        dualVolume = rhs.dualVolume;
        return (*this);
    }

    int Get ( int i ) const {return v[i]; }
    unsigned int operator[] (int i)const {
        return v[i];
    }
    unsigned int& operator[] (int i) {
        return v[i];
    }
    friend inline std::ostream& operator<<(std::ostream& os, const Simplex<NT,Dim> & simplex)
    {
        os << "(" << simplex.Index() << ": ";
        if(simplex.IsNegative())
        {
            os << "-";
        }
        else
        {
            os << "+";
        }
        os << simplex.Volume();

        os << simplex.getIndexSet() << ")";
        return os;
    }


protected:
    int sign=1;
    size_t index=-1;
    Vector center=Vector::Zero();
    Scalar volume=0;
    Scalar dualVolume=0;
    mtao::IndexSet<Dim+1> v;



    /*
     * OPERATORS
     *
     */
public:

    inline bool operator==(  const Simplex<NT,Dim>& other ) const
    {
        return this->getIndexSet() == other.getIndexSet();
    }

    inline bool operator!=(  const Simplex<NT,Dim>& other ) const
    {
        return (!(*this == other));
    }

    inline bool operator<(  const Simplex<NT,Dim>& other ) const
    {
        return this->getIndexSet() < other.getIndexSet();
    }

    inline bool operator>(  const Simplex<NT,Dim>& other ) const
    {
        return (other < *this);
    }

    inline bool operator<=(  const Simplex<NT,Dim>& other ) const
    {
        return (!(other < this));
    }

    inline bool operator>=(  const Simplex<NT,Dim>& other ) const
    {
        return (!(*this < other));
    }

    template <int Dim2>
    inline bool isSameSign( const Simplex<NT,Dim2>& other ) const
    {
        return (this->IsNegative() == other.IsNegative());
    }


};
template<typename NT, int N>
inline Simplex<NT,N> operator-( const Simplex<NT,N>& simplex )
{
    Simplex<NT,N> negated( simplex );
    negated.Negate();
    return negated;
}



template <typename NT>
inline Simplex<NT,0> MakeSimplex( int v0 )
{
    int v[1] = {v0};
    return Simplex<NT,0>( v );
}

template <typename NT>
inline Simplex<NT,1> MakeSimplex( int v0, int v1 )
{
    int v[2] = {v0, v1};
    return Simplex<NT,1>( v );
}

template <typename NT>
inline Simplex<NT,2> MakeSimplex( int v0, int v1, int v2 )
{
    int v[3] = {v0, v1, v2};
    return Simplex<NT,2>( v );
}

template <typename NT>
inline Simplex<NT,3> MakeSimplex( int v0, int v1, int v2, int v3 )
{
    int v[4] = {v0, v1, v2, v3};
    return Simplex<NT,3>( v, false );
}

#endif
