#include "../include/trianglemesh.h"
#include "../include/advection.hpp"
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
    return std::move(normals);
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
    return std::move(normals);
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
    return std::move(normals);
}
auto NormalTriangleMesh:: meanCurvatureNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));


    for(auto&& s: simplices<2>()) {
        for(typename SparseMatrix::InnerIterator it(b(), s.Index()); it; ++it) {
            auto& e = simplex<1>(it.row());
            uint ind = s.oppositeIndex(e);
            Vector v0 = vertex(e[0]) - vertex(ind);
            Vector v1 = vertex(e[1]) - vertex(ind);
            Scalar cot = v0.dot(v1)/v0.cross(v1).norm();//If the mesh is voronoi all angles should be less than \pi/2 right?
            Vector dir = cot * vertex(e[1]) - vertex(e[0]);
            normals[e[0]] -= dir;
            normals[e[1]] -= dir;
        }
    }
    int i=0;
    for(auto&& n: normals) {
        EIGEN_DEBUG_VAR(vertex(i).normalized().dot(n))
        n.normalize();
    }

    return std::move(normals);
}
auto NormalTriangleMesh:: sphereNormal() -> std::vector<Vector> {
    std::vector<Vector> normals(vertices().size(), Vector(0,0,0));
    return std::move(normals);
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




void NormalTriangleMesh::semilagrangianAdvection(VelocityFormType & form, Scalar dt){
    VelocityFormType oldform = form;
    std::cout << oldform.expr.rows() << " " << form.expr.rows() << std::endl;

    for(auto&& s2: constSimplices<2>()) {
        typename SparseMatrix::InnerIterator it2(b(), s2.Index());
        int col = 0;
        if(it2.value() * form(it2.row()) < 0 ) {//==1 => inward vector, which we dont want because we're back advecting
            ++it2;
            col = 1;
        }
        if(it2) {
            auto&& s1 = simplex<1>(it2.row());
            Particle<NormalTriangleMesh> p(*this,s1.Center(), &s2);
            p.advectInPlace(oldform,-dt);
            Vector newvel = getVelocity(p.p(),s2,oldform);
            form.expr(it2.row()) = it2.value() * newvel.dot(whitneyBasis(s2).col(col));
        } else {
            std::cout << "WTF i just intersected with a boundary" << std::endl;
        }


    }
    std::cout << (form - oldform).expr.norm() << std::endl;
}

void NormalTriangleMesh::advect(VelocityFormType & form, Scalar dt, AdvectionType type){
    switch(type) {
    case Semilagrangian_Advection:
        semilagrangianAdvection(form,dt);
    }
}
