#ifndef DEC_GEOMETRY_H
#define DEC_GEOMETRY_H
#include "dec.hpp"
#include <vector>
#include <limits>
#include <Eigen/Geometry>
namespace mtao{

constexpr int cefactorial(int n)
{
    return n > 0 ? n * cefactorial(n-1):1;
}

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
void normalizeToBBoxInPlace(Vector & v, const Eigen::AlignedBox<typename Vector::Scalar, Vector::RowsAtCompileTime> & bbox) {
    typedef typename Vector::Scalar Scalar;
    Vector mid = bbox.center();
    Scalar range = (bbox.max()-bbox.min()).maxCoeff();
    v.noalias() = (v-mid)/range;
}



template <typename Vector>
void unnormalizeToBBoxInPlace(Vector & v, const Eigen::AlignedBox<typename Vector::Scalar, Vector::RowsAtCompileTime> & bbox) {
    typedef typename Vector::Scalar Scalar;
    Vector mid = bbox.center();
    Scalar range = (bbox.max()-bbox.min()).maxCoeff();
    v.noalias() = (v+mid)*range;
}








template <typename Vector>
const std::vector<Vector> normalizeToBBox(const std::vector<Vector> & vertices, const Eigen::AlignedBox<typename Vector::Scalar, Vector::RowsAtCompileTime> & bbox) {
    std::vector<Vector> ret = vertices;
    normalizeToBBoxInPlace(ret,bbox);
    return ret;
}
template <typename Vector>
Vector normalizeToBBox(const Vector & vertices, const Eigen::AlignedBox<typename Vector::Scalar, Vector::RowsAtCompileTime> & bbox) {
    Vector ret = vertices;
    normalizeToBBoxInPlace(ret,bbox);
    return ret;
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



/*
template <typename Scalar, bool Signed>
Scalar volume(const Eigen::Matrix<Scalar,Eigen::Dynamic,Eigen::Dynamic> & m) {
    if(Signed) {
        return m.determinant()/cefactorial(Dim);
    } else {
        return std::sqrt((m.transpose()*m).determinant())/cefactorial(Dim);
    }
}
*/

template <typename Scalar,int EmbeddedDim, int Dim>
Scalar signedVolume(const Eigen::Matrix<Scalar,EmbeddedDim,Dim> & m) {
    return m.determinant()/cefactorial(Dim);
}

template <typename Scalar,int EmbeddedDim, int Dim>
Scalar unsignedVolume(const Eigen::Matrix<Scalar,EmbeddedDim,Dim> & m) {
    return std::sqrt((m.transpose()*m).determinant())/cefactorial(Dim);
}
template <typename Scalar,int EmbeddedDim, int Dim, bool Signed = true>
Scalar volume(const Eigen::Matrix<Scalar,EmbeddedDim,Dim> & m) {
    if(Dim == EmbeddedDim && Signed) {
        return signedVolume(m);
    } else {
        return unsignedVolume(m);
    }
}
















template <typename Matrix>
void gramSchmidt(Matrix & basis) {
    for(int i=0; i < basis.cols(); ++i) {
        basis.col(i).normalize();
        basis.rightCols(basis.cols()-i-1) -=  basis.col(i) *  (basis.col(i).transpose() * basis.rightCols(basis.cols()-i-1));
    }
    basis.col(basis.cols()-1).normalize();
}

template <typename Matrix>
auto normal(Matrix & basis) -> Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1>{
    typedef Eigen::Matrix<typename Matrix::Scalar, Matrix::RowsAtCompileTime, 1> Vector;
    gramSchmidt(basis);

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
