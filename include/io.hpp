#ifndef MESH_UTIL_H
#define MESH_UTIL_H
#include <boost/tokenizer.hpp>
#include <boost/algorithm/string.hpp>
#include "simplicialComplex.hpp"
#include <fstream>


template <typename Scalar>
SimplicialComplex<NumericalTraits<Scalar,3>,2> * readOBJtoSimplicialComplex(std::istream & is)
{

    typedef SimplicialComplex<NumericalTraits<Scalar,3>, 2> TriangleMesh;
    typedef Eigen::Matrix<Scalar,3,1> Vector3;
    std::vector<Vector3> verts;
    std::vector<mtao::IndexSet<3> > indexsets;
    //temporary variables
    std::string line;
    std::string s;
    Vector3 x;
    mtao::IndexSet<3> f;
#ifdef CHECK_FOR_ZERO_OBJ
    unsigned int minVertIndex=100000;
#endif

    //expect lines to only be up to this long
    line.reserve(64);
    while( is.good() )
    {
        getline(is,line);
        if(line.length()>0)
        {
            boost::trim(line);
            std::stringstream ss(line);
            std::string type;
            ss >> type;

            if(type.compare("v") == 0)
            {
                //vertex
                ss >> x(0) >> x(1) >> x(2);
                if(!ss.eof())
                {
                    Scalar w;
                    ss >> w;
                    x/=w;
                }
                verts.push_back(x);
            }
            else if(type[0]=='f')
            {
                //Face (assumed to be triangular or quad)
                for( int i=0; i<3; ++i)
                {

                    ss >> s;
                    boost::tokenizer<> tok(s);

                    boost::tokenizer< >::iterator tok_it = tok.begin();
                    if(tok_it!=tok.end())
                    {//vertex
                        f[i] = atoi(tok_it->c_str())-1;
#ifdef CHECK_FOR_ZERO_OBJ
                        minVertIndex=(minVertIndex>f[i]+1)?f[i]+1:minVertIndex;
#endif

                    }
                }
                indexsets.push_back(f);
                if(!ss.eof())//maybe this is a quad, if it's more than a quad I won't deal with it (and just fail)
                {

                    ss >> s;
                    boost::tokenizer<> tok(s);

                    boost::tokenizer< >::iterator tok_it = tok.begin();
                    if(tok_it!=tok.end())
                    {//vertex
                        //swap first and last to keep ccw
                        f[0]^=f[2];
                        f[2]^=f[0];
                        f[0]^=f[2];
                        f[1] = atoi(tok_it->c_str())-1;
#ifdef CHECK_FOR_ZERO_OBJ
                        minVertIndex=(minVertIndex>f[1]+1)?f[1]+1:minVertIndex;
#endif

                    }
                    indexsets.push_back(f);


                }

            }
            else
            {
                //Comment or something we don't know how to parse

            }

        }

    }

#ifdef CHECK_FOR_ZERO_OBJ
    if(minVertIndex==0)
    {
        std::cout << "This wasn\'t a real OBJ file, vertices started at 0\n";
        for(auto && indexset: indexsets)
        {
            indexset[0]+=1;
            indexset[1]+=1;
            indexset[2]+=1;
        }
    }
#endif
    return new TriangleMesh(indexsets, verts);

}


template <typename Scalar>
SimplicialComplex<NumericalTraits<Scalar,3>,2> * readOBJtoSimplicialComplex(const std::string & str)
{
    std::ifstream is(str.c_str());
    if(!is.is_open()){return 0;}
    return readOBJtoSimplicialComplex<Scalar>(is);
}

template <typename SimplicialComplex>
void writeSimplicialComplextoStream(const SimplicialComplex & sc, std::ostream & os)
{
    static_assert(SimplicialComplex::Dim > 0, "No point in pulling from a 0 complex...");
    const unsigned int Dim = SimplicialComplex::NumTraits::Dim;

    for(auto && v: sc.constVertices())
    {
        os << "v "<< v.transpose() << std::endl;
    }
    for(auto && s: sc.constSimplices())
    {
        os << "f ";
        if(s.IsNegative()) {
            os << s[1] << " " << s[0] << " ";
        }
        else
        {
            os << s[0] << " " << s[1] << " ";
        }
        for(int i=2; i < Dim; ++i)
        {
            os << s[i] << " ";
        }
        os << std::endl;
    }
}

template <typename SimplicialComplex>
void writeOBJfromSimplicialComplex(const SimplicialComplex & sc,const std::string & str)
{
    std::cout << "Writing SC to " << str << std::endl;
    std::ofstream os(str.c_str());
    writeSimplicialComplextoStream(sc,os);
}


#endif
