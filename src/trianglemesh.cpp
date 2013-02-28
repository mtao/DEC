#include "../include/trianglemesh.h"
#include <cmath>
auto NormalTriangleMesh::equalNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));



    for(auto&& s: simplices<2>()) {
        for(int i=0; i < 3; ++i)
        {
        normals[s[i]] +=  m_face_normals[s.Index()];
        }
    }
    for(auto&& n: normals) {
        n.normalize();
    }
    return normals;
}

auto NormalTriangleMesh:: areaNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));
    for(auto&& s: simplices<2>()) {
        for(int i=0; i < 3; ++i)
        {
        normals[s[i]] += s.Volume() * m_face_normals[s.Index()];
        }
    }
    for(auto&& n: normals) {
        n.normalize();
    }
    return normals;
}
auto NormalTriangleMesh:: angleNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));
    for(auto&& s: simplices<2>()) {
        for(int i=0; i < 3; ++i)
        {
            int other[2];
            other[0] = (i+1)%3;
            other[1] = (i+2)%3;
            auto&& origin = vertex(s[i]);
            auto&& left = (vertex(s[other[0]]) - origin).normalized();
            auto&& right = (vertex(s[other[1]]) - origin).normalized();
            Scalar factor = std::acos((left.dot(right)));
        normals[s[i]] += factor * m_face_normals[s.Index()];
        }
    }
    for(auto&& n: normals) {
        n.normalize();
    }
    return normals;
}
auto NormalTriangleMesh:: meanCurvatureNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));
    return normals;
}
auto NormalTriangleMesh:: sphereNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));
    return normals;
}
auto NormalTriangleMesh:: getNormals(NormalType type) -> std::vector<Vector> {
    switch(type) {
    case Equal_Normal:
        return equalNormal();
        break;
    case Area_Normal:
        return areaNormal();
        break;
    case Angle_Normal:
        return angleNormal();
        break;
    case Mean_Curvature_Normal:
        return meanCurvatureNormal();
        break;
    case Sphere_Normal:
        return sphereNormal();
        break;
    }
}

void NormalTriangleMesh::init() {
    computeNormals();
}

void NormalTriangleMesh::computeNormals() {
    m_face_normals.resize(simplices().size());
    for(auto&& s: simplices<2>()) {
        Eigen::Matrix3f mat = vertices(s);
        mat.leftCols(2) = mat.leftCols(2) - mat.col(2).rowwise().replicate(2);
        m_face_normals[s.Index()] =
                (
                (s.isNegative()?-1:1) *
                 mat.col(0).cross(mat.col(1)).normalized());
    }
}
