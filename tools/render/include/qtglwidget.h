#ifndef GLWIDGET_H
#define GLWIDGET_H
#include "glutil.h"
#include <QGLWidget>
#include <glm/glm.hpp>
#include <memory>
#include <QTimer>
#include <QTime>
#include <QVector3D>

struct MeshPackage {
    std::shared_ptr<VertexBufferObject> vertices;
    std::shared_ptr<VertexIndexObject> indices;
    std::shared_ptr<VertexBufferObject> facevertices;
    std::shared_ptr<VertexIndexObject> faceindices;
    std::shared_ptr<VertexBufferObject> edgevertices;
    std::shared_ptr<VertexIndexObject> edgeindices;
};


enum RenderType {RT_FACE=4, RT_VERT=1, RT_EDGE=2, RT_NONE=0};

struct FormPackage{
    QString title;
    RenderType type;
    std::shared_ptr<VertexBufferObject> data;
};

class MainWindow;
class GLWidget: public QGLWidget
{
    Q_OBJECT
public:
    friend class MainWindow;
    GLWidget(QWidget * parent);
protected:
    void paintGL();
    void resizeGL(int w, int h);
    virtual void initializeGL();
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * wheelEvent);
    virtual void keyPressEvent(QKeyEvent *event);
private:
    void initShader(ShaderProgram & program, const QString & geotype);
    GLuint compileShader(GLenum shaderType, const QString & fileName);
    std::unique_ptr<ShaderProgram> & shaderSelector(RenderType type);
    void render(RenderType type);
    //    void renderForm(const FormPackage & form);
    QTimer* m_timer = 0;
    int m_renderType = RT_FACE;
    bool m_doRender = false;
    glm::mat4 mat_mvp, mat_mv, mat_m, mat_v, mat_p;
    float m_aspectRatio = 1;
    std::unique_ptr<ShaderProgram> m_shader;
    std::unique_ptr<ShaderProgram> m_vertshader;
    std::unique_ptr<ShaderProgram> m_faceshader;
    std::unique_ptr<ShaderProgram> m_edgeshader;
    MeshPackage m_meshpackage;
    std::map<QString, FormPackage> m_formpackages;

    QTime m_time;
    int m_lastTime = 0;
    int m_mouseEventTime;
    float m_distance = 1.4f;
    QPointF lastPos;
    QVector3D m_rotation;
    QVector3D m_angularMomentum;
    QVector3D m_accumulatedMomentum;
public slots:
    void disableRendering(){m_doRender =false;}
    void enableRendering(){m_doRender = true;}
    void recieveMesh(const MeshPackage & package);
    void recieveForm(const FormPackage & package);
    void enableForm(const QString & formname);
};
#endif
