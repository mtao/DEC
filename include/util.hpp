#ifndef DEC_UTILS_H
#define DEC_UTILS_H
#include "dec.hpp"
#include <vector>
#include <limits>
#include "geometry.hpp"
namespace mtao{



template <typename Complex>
void normalInPlace(const Complex & sc, typename Complex::Vector & n, const typename Complex::template TraitsContainer<Complex::Dim>::simplextype & simplex) {
    if(Complex::EmbeddedDim == 3 && decltype(simplex)::Dim == 2) {
        n = (simplex.isNegative()?1:-1)*(sc.vertex(simplex[2]) - sc.vertex(simplex[0])).cross(sc.vertex(simplex[2]) - sc.vertex(simplex[0])).normalize();
    }
}
template <typename Complex>
void projectToSimplexInPlace(const Complex & sc, typename Complex::Vector & vec, const typename Complex::template TraitsContainer<Complex::Dim>::simplextype & simplex) {

    typedef typename Complex::Vector Vector;
    auto&& origin = sc.vertex(simplex[Complex::Dim]);
    Vector ortho = vec - origin;
    Eigen::Matrix<typename Vector::Scalar, Vector::RowsAtCompileTime, Complex::Dim> m;
    for(int i=0; i < Complex::Dim; ++i) {
        m.col(i) = sc.vertex(simplex[i]) - origin;
    }
    Vector n = normal(m);
    ortho = ortho.dot(n) * n;
    //ortho.normalize();

    vec -= ortho;//.dot(vec) * ortho;

}
template <typename Complex>
typename Complex::Vector projectToSimplex(const Complex & sc, const typename Complex::Vector & vec, const typename Complex::template TraitsContainer<Complex::Dim>::simplextype & simplex) {
    typename Complex::Vector v = vec;
    projectToSimplexInPlace(sc,v,simplex);
    return v;


}


};

#endif
