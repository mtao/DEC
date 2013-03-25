#include <iostream>
#include "simplicialComplex.hpp"
#include "io.hpp"
#include "dec.hpp"

class DECTest: public DEC<TriangleMesh>
{
    public:
    DECTest(const TriangleMesh & sc): DEC<TriangleMesh>(sc)
    {
        std::cout << "Derivative: " << std::endl;
        std::cout << d<0>().constData() << std::endl;
        std::cout << d<1>().constData() << std::endl;
        std::cout << "Hodges: " << std::endl;
        std::cout << h<0>().constData().diagonal().transpose() << std::endl;
        std::cout << h<1>().constData().diagonal().transpose() << std::endl;
        std::cout << h<2>().constData().diagonal().transpose() << std::endl;

        std::cout << "Derivative product: " << std::endl;
        std::cout << (d<1>()*d<0>()).constData() << std::endl;
        std::cout << "Derivatie composition: " << std::endl;
        std::cout << d(d<0>()).constData() << std::endl;
        std::cout << "Laplace downward" << std::endl;
        std::cout << h(d(h(d<0>()))).constData() << std::endl;
        std::cout << "Laplace1" << std::endl;
        std::cout << (h(d(h(d<1>()))) + d(h(d(h<1>())))).constData() << std::endl;
        std::cout << "Laplace2" << std::endl;
        std::cout << ((d(h<2>()))).constData() << std::endl;
        std::cout << d(h(d(h<2>()))).constData() << std::endl;
     auto form = genForm<FormType::Primal,1>();
    auto form2 = genForm<FormType::Primal,1>();
    form2.data().setOnes();
    decltype(genForm<FormType::Primal,1>()) form3(form2 + form2);
    form3 = form3 + form2;
    std::cout << form.constData().transpose() << std::endl;
    std::cout << form2.constData().transpose() << std::endl;
    std::cout << form3.constData().transpose() << std::endl;
    std::cout << d(form3).constData() << std::endl;
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
