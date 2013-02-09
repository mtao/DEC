#ifndef DEC_UTILS_H
#define DEC_UTILS_H
#include "dec.hpp"
#include <vector>
#include <limits>
namespace mtao{

template <typename Vector>
void normalizeInPlace(std::vector<Vector> & vertices) {
    typedef typename Vector::Scalar Scalar;
    Vector min = Vector::Constant(std::numeric_limits<Scalar>::infinity());
    Vector max = Vector::Constant(-std::numeric_limits<Scalar>::infinity());
    for(auto&& vert: vertices) {
        min = min.cwiseMin(vert);
        max = max.cwiseMax(vert);
    }
    Vector mid = (max+min)/2;
    Scalar range = (max-min).maxCoeff();
    for(auto&& v: vertices) {
        v.noalias() = (v-mid)/range;
    }
}
template <typename Vector>
std::vector<Vector> normalize(const std::vector<Vector> & vertices) {
    std::vector<Vector> ret = vertices;
    normalizeInPlace(ret);
    return ret;
}

};

#endif
