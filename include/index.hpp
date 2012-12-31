#ifndef INDEX_UTIL_H
#define INDEX_UTIL_H
#include <algorithm>

namespace mtao
{
template <unsigned int N, typename Index=unsigned int>
struct IndexSet{
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
    Index m_data[N];
    Index & operator[](int i)
    {
        return m_data[i];
    }
    Index Get ( int i ) const {return m_data[i];}
    const Index & operator[](int i) const
    {
        return m_data[i];
    }
    IndexSet()
    {
        std::fill(m_data, m_data+N, Index(0));
    }

    IndexSet(const IndexSet & other)
    {
        std::copy(other.m_data, other.m_data+N, m_data);
    }

    IndexSet(const Index ind[N])
    {
        std::copy(ind,ind+N, m_data);
    }
    IndexSet(std::initializer_list<Index> ind)
    {
        std::copy(ind.begin(),ind.end(), m_data);
    }


    const static unsigned int dim = N;

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
                return true;
            else if( this->Get(i) > other.Get(i) )
                return false;
        }

        return (this->Get(N) < other.Get(N));
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
