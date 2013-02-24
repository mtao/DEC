#ifndef DEC_GEOMETRY_H
#define DEC_GEOMETRY_H
#include "dec.hpp"
#include <vector>
#include <limits>
#include <Eigen/Geometry>
namespace mtao{

template <typename Vector>
void getBoundingBoxInPlace(const std::vector<Vector> & vertices, Eigen::AlignedBox<typename Vector::Scalar,Vector::RowsAtCompileTime> & bbox)  {
    typedef typename Vector::Scalar Scalar;
    auto&& min = bbox.min();
    auto&& max = bbox.max();
    min = Vector::Constant(std::numeric_limits<Scalar>::infinity());
    max = Vector::Constant(-std::numeric_limits<Scalar>::infinity());
    for(auto&& vert: vertices) {
        min = min.cwiseMin(vert);
        max = max.cwiseMax(vert);
    }
}

template <typename Vector>
Eigen::AlignedBox<typename Vector::Scalar,Vector::RowsAtCompileTime> getBoundingBox(const std::vector<Vector> & vertices) {
    Eigen::AlignedBox<typename Vector::Scalar,Vector::RowsAtCompileTime> bbox;
    getBoundingBoxInPlace(vertices,bbox);
    return bbox;
}

template <typename Vector>
void normalizeToBBoxInPlace(std::vector<Vector> & vertices, const Eigen::AlignedBox<typename Vector::Scalar, Vector::RowsAtCompileTime> & bbox) {
    typedef typename Vector::Scalar Scalar;
    Vector mid = bbox.center();
    Scalar range = (bbox.max()-bbox.min()).maxCoeff();
    for(auto&& v: vertices) {
        v.noalias() = (v-mid)/range;
    }
}

template <typename Vector>
void normalizeInPlace(std::vector<Vector> & vertices) {
    auto bbox = getBoundingBox(vertices);
    normalizeToBBoxInPlace(vertices,bbox);
}

template <typename Vector>
std::vector<Vector> normalize(const std::vector<Vector> & vertices) {
    std::vector<Vector> ret = vertices;
    normalizeInPlace(ret);
    return ret;
}


template <typename Matrix>
void orthogonalize(Matrix & basis) {
    for(int i=0; i < basis.cols(); ++i) {
        basis.col(i).normalize();
        for(int j=i+1; j < basis.cols(); ++j) {
            basis.col(j) -= basis.col(i).dot(basis.col(j)) * basis.col(i);
        }
    }
}

template <typename Matrix>
auto normal(Matrix & basis) -> Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1>{
    typedef Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1> Vector;
    typedef typename Vector::Scalar Scalar;
    orthogonalize(basis);

    Vector v;
    while(true){
        v = Vector::Random();
        v -= basis * (basis.transpose() * v);
        /*
        for(int i=0; i < basis.cols(); ++i) {
            v -= basis.col(i).dot(v) * v.col(i);
        }
        */
        if(v.norm() > 0.0001) {
            return v.normalized();
        }
    }


}

template <typename Matrix>
auto constNormal(const Matrix & mat) -> Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1>{
    typedef Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1> Vector;
    typedef typename Vector::Scalar Scalar;

        Eigen::Matrix<Scalar, Matrix::RowsAtCompileTime, Matrix::ColsAtCompileTime> basis(mat);
        return normal(basis);


}


};

#endif
