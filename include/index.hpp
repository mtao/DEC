#ifndef INDEX_UTIL_H
#define INDEX_UTIL_H
#include <algorithm>
#include <array>

namespace mtao
{
template <unsigned int N, typename Index=unsigned int>
struct IndexSet: public std::array<Index, N>{
    friend inline std::ostream& operator<<(std::ostream& os, const mtao::IndexSet<N,Index> & is)
    {
        os << "[" << is.Get(0);
        for(unsigned int i=1; i < N; ++i)
        {
            os << "," << is.Get(i);
        }
        os << "]";
        return os;
    }
    Index Get ( int i ) const {return this->at(i);}
    IndexSet()
    {
        this->fill(0);
    }

    IndexSet(const IndexSet & other)
    {
        std::copy(other.cbegin(), other.cend(), this->begin());
    }

    IndexSet(const Index ind[N])
    {
        std::copy(ind,ind+N, this->begin());
    }
    IndexSet(std::initializer_list<Index> ind)
    {
        std::copy(ind.begin(),ind.end(), this->begin());
    }



public:
    inline bool operator==(const IndexSet<N,Index>& other ) const
    {
        for( unsigned int i = 0; i < N; ++i )
            if( this->Get(i) != other.Get(i) )
                return false;
        return true;
    }

    inline bool operator!=( const IndexSet<N,Index>& other ) const
    {
        return (!(*this == other));
    }

    inline bool operator<( const IndexSet<N,Index>& other ) const
    {
        for( unsigned int i = 0; i < N; ++i )
        {
            if( this->Get(i) < other.Get(i) )
            {
                return true;
            }
            else if( this->Get(i) > other.Get(i) )
            {
                return false;
            }
        }

        return false;
    }

    inline bool operator>( const IndexSet<N,Index>& other ) const
    {
        return (other < *this);
    }

    inline bool operator<=( const IndexSet<N,Index>& other ) const
    {
        return (!(other < *this));
    }

    inline bool operator>=( const IndexSet<N,Index>& other ) const
    {
        return (!(*this < other));
    }
};
}

#endif
