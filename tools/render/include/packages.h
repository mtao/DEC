
#ifndef PACKAGES_H
#define PACKAGES_H
#include <QString>
#include <memory>
#include "glutil.h"
#include <Eigen/Dense>
struct MeshPackage {
    const std::vector<Eigen::Vector3f> vertices;
    const std::vector<unsigned int> triangle_indices;
    const std::vector<unsigned int> edge_indices;
    const std::vector<Eigen::Vector3f> dual_vertices;
    const std::vector<unsigned int> dual_edge_indices;
    const std::vector<unsigned int> dual_vertex_indices;

    /*
    std::shared_ptr<VertexBufferObject> vertices;
    std::shared_ptr<VertexIndexObject> indices;
    std::shared_ptr<VertexBufferObject> facevertices;
    std::shared_ptr<VertexIndexObject> faceindices;
    std::shared_ptr<VertexBufferObject> edgevertices;
    std::shared_ptr<VertexIndexObject> edgeindices;
    */
    /*
    std::shared_ptr<VertexBufferObject> dualfacevertices;
    std::shared_ptr<VertexIndexObject> dualfaceindices;
    std::shared_ptr<VertexBufferObject> dualedgevertices;
    std::shared_ptr<VertexIndexObject> dualedgeindices;
    */
};


enum RenderType {RT_FACE=4, RT_VERT=1, RT_EDGE=2, RT_NONE=0};

struct FormPackage{
    QString title;
    RenderType type;
    std::shared_ptr<VertexBufferObject> data;
};



template <typename Form>
constexpr RenderType formToRendertype(const Form &) {
    return RenderType((Form::Traits::NOut == 0)
                      *RT_VERT+
                      (Form::Traits::NOut == 1)
                      * RT_EDGE+
                      (Form::Traits::NOut == 2)
                      * RT_FACE);


}
template <typename Form>
FormPackage makeFormPackage(const QString & name, const Form & form)
{
    auto&& tmparr = formToRenderable(form);
    return  {name,formToRendertype(form), std::make_shared<VertexBufferObject>(
                tmparr.data(),
                (Form::Traits::NOut+1)*tmparr.size(),
                GL_STATIC_DRAW,1)};
}
#endif
