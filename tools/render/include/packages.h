#ifndef PACKAGES_H
#define PACKAGES_H
#include <memory>
#include <QString>
#include "glutil.h"
#include <Eigen/Dense>
#include "../../../include/dec.hpp"

struct MeshBuffers {
    std::shared_ptr<VertexBufferObject> vertices;
    std::shared_ptr<VertexIndexObject> indices;
    std::shared_ptr<VertexBufferObject> facevertices;
    std::shared_ptr<VertexIndexObject> faceindices;
    std::shared_ptr<VertexBufferObject> edgevertices;
    std::shared_ptr<VertexIndexObject> edgeindices;
    std::shared_ptr<VertexBufferObject> dual_vertices;
    std::shared_ptr<VertexIndexObject> dual_indices;
    std::shared_ptr<VertexBufferObject> dual_facevertices;
    std::shared_ptr<VertexIndexObject> dual_faceindices;
    std::shared_ptr<VertexBufferObject> dual_edgevertices;
    std::shared_ptr<VertexIndexObject> dual_edgeindices;
};
struct MeshPackage {
    std::vector<Eigen::Vector3f> vertices;
    std::vector<unsigned int> indices;
    std::vector<Eigen::Vector3f> facevertices;
    std::vector<unsigned int> faceindices;
    std::vector<Eigen::Vector3f> edgevertices;
    std::vector<unsigned int> edgeindices;
    std::vector<Eigen::Vector3f> dual_vertices;
    std::vector<unsigned int> dual_indices;
    std::vector<Eigen::Vector3f> dual_facevertices;
    std::vector<unsigned int> dual_faceindices;
    std::vector<Eigen::Vector3f> dual_edgevertices;
    std::vector<unsigned int> dual_edgeindices;
};


enum RenderType {RT_FACE=4, RT_VERT=1, RT_EDGE=2, RT_NONE=0, RT_DUAL = 8};

struct FormPackage{
    QString title;
    RenderType type;
    std::shared_ptr<VertexBufferObject> data;
};

namespace mtao{
    template <typename Form>
        constexpr RenderType formToRendertype(const Form &) {
            return RenderType(
                        (Form::Traits::TypeOut== DUAL_FORM)
                *RT_DUAL+
                        (Form::Traits::NOut == 0)
                *RT_VERT+
            (Form::Traits::NOut == 1)
                * RT_EDGE+
            (Form::Traits::NOut == 2)
                * RT_FACE);


        }
    template <typename Form>
        FormPackage makeFormPackage(const QString & name, const Form & form, const std::vector<unsigned int> & indices = std::vector<unsigned int>())
        {
            auto&& tmparr = formToRenderable(form, indices);
            return  {name,formToRendertype(form), std::make_shared<VertexBufferObject>(
                    tmparr.data(),
                    (Form::Traits::NOut+1)*tmparr.size(),
                    GL_STATIC_DRAW,1)};
        }
};
#endif
