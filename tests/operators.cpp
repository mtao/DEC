#include <iostream>
#include "simplicialComplex.hpp"
#include "io.hpp"
#include "dec.hpp"

class DECTest: public DEC<TriangleMesh>
{
    public:
    DECTest(const SimplicialComplex & sc): DEC<TriangleMesh>(sc)
    {
        std::cout << d<1>().constData() << std::endl;
    Form<typename SimplicialComplex::NumTraits::DynamicVector,PRIMAL_FORM,1> form = genForm<PRIMAL_FORM,1>();
    auto form2 = genForm<PRIMAL_FORM,1>();
    form2.data().setOnes();
    Form<typename SimplicialComplex::NumTraits::DynamicVector,PRIMAL_FORM,1> form3(form2 + form2);
    std::cout << form.constData().transpose() << std::endl;
    std::cout << form2.constData().transpose() << std::endl;
    std::cout << form3.constData().transpose() << std::endl;
    std::cout << "Done! " << std::endl;
    }

};

int main()
{
    typedef Eigen::Vector3d V;
    std::vector<V> vecs;
    std::vector<mtao::IndexSet<3> > tris;
    for(int i=0; i < 10; ++i)
    {
        vecs.push_back(V(i,i,i));
    }
    vecs[0] = V(0,0,0);
    vecs[3] = V(1,1,0);
    vecs[5] = V(1,0,0);
    vecs[6] = V(0,1,0);
    tris.push_back(mtao::IndexSet<3>({0,3,5}));
    tris.push_back(mtao::IndexSet<3>({3,0,6}));

    TriangleMesh sc(tris,vecs);
    DECTest dec(sc);


    return 0;
}