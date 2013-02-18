#ifndef GLUTIL_H
#define GLUTIL_H

#include <GL/glew.h>
#include <vector>
#include <map>
#include <memory>

/*
 Standard use:
initializeGL()
{
ShaderProgram p;
VertexBufferObject vbo1;
VertexBufferObject vbo2;
VertexBufferObject vbo3;
p.addAttribute("vertex",&vbo1);
p.addAttribute("normal",&vbo2);
p.addAttribute("color",&vbo3);
VertexIndexObject vio;
}
paintGL()
{
program.bind();
vio.render();
program.release();
}

  */


struct VertexBufferObject
{
    /*
      Capabilities:
      create vbo
      load vbo
      */
    VertexBufferObject(): owner(false) {}
    VertexBufferObject(GLvoid *data, GLsizei size, GLenum usage=GL_STATIC_DRAW,
                       int tupleSize=3, GLenum type=GL_FLOAT, GLboolean normalized=GL_FALSE,
                       GLsizei stride=0, GLvoid * offset=0);
    void bind(GLint attributeId);
    ~VertexBufferObject();
    GLuint vbo;
    int tupleSize;//Values stored redundantly between constructor but doesn't matter
    int size;
    GLenum type;
    GLboolean normalized;
    GLsizei stride;
    GLvoid * offset;
    bool owner=true;
};
struct VertexIndexObject
{
    /*
      Capabilities:
      create ibo
      render using ibo
      */
    VertexIndexObject(): owner(false) {}
    VertexIndexObject(GLvoid *data, GLsizei size, GLenum usage,
                      GLenum mode=GL_TRIANGLES, GLenum type=GL_UNSIGNED_INT, GLvoid * offset=0);
    ~VertexIndexObject();
    void render();
    GLuint ibo;
    GLenum mode;//Once again values stored redundantly between constructor but doesn't matter
    GLsizei size;
    GLenum type;
    GLvoid * offset;
    bool owner=true;
};

struct ShaderProgram
{
    /*
      Capabilities:
      Create shader program
      bind shader program (and perhaps all of the attributes it has bound)
      release shader program (and perhaps all of the attributes it has bound)

      optionally:
      check for attributes
      check for redundant attributes
      */
    ShaderProgram();
    ShaderProgram(GLuint programId);
    void init();
    void addAttribute(const char * name, std::shared_ptr<VertexBufferObject> vbo);
    GLint getAttribLocation(const char * name);
    std::shared_ptr<VertexBufferObject> getAttribute(const char * name);
    void bind(bool bindAttributes=true);
    void release();
    GLuint programId;
    bool attributesBound=false;
    //std::map<GLint, std::weak_ptr<VertexBufferObject> > attributes;
    std::map<GLint, std::shared_ptr<VertexBufferObject> > attributes;
    //std::vector< VertexBufferObject> attributes;
    bool owner=true;
};


struct ShaderCollection
{
    /*
      Capabilities:
      Holds a shader program, its attributes, and an ibo
      */
    ShaderCollection() {}
    void setProgram(ShaderProgram & program_) {program=program_; program.owner=false;}
    void setVertexIndexObject(VertexIndexObject & vio_) {vio=vio_; vio.owner=false;}
    void addVertexBufferObject(VertexBufferObject & vbo) {vbos.push_back(vbo); vbos.back().owner=false;}
    void render(){vio.render();}
    ShaderProgram program;
    VertexIndexObject vio;
    std::vector<VertexBufferObject> vbos;

};

#endif // GLUTIL_H
