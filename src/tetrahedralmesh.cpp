#include "../include/tetrahedralmesh.h"



auto TetrahedralMesh::getGradient(int triidx,int tetidx) const -> Vector {
    auto&& d = this->template d<3>();
    assert(triidx >= 0 && triidx < d.cols());
    assert(tetidx >= 0 && tetidx < d.rows());
    float sgn = d.coeff(tetidx,triidx);
    auto&& tet = cell(tetidx);
    auto&& tri = face(triidx);
    
    const Vector& v0 = vertex(tri[0]);
    const Vector& v1 = vertex(tri[1]);
    const Vector& v2 = vertex(tri[2]);

    auto&& a = v1 - v0;
    auto&& b = v2 - v0;

    Vector normal = sgn * a.cross(b);
}
