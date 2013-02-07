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
};

class GLWidget: public QGLWidget
{
    Q_OBJECT
public:
    GLWidget(QWidget * parent);
protected:
    void paintGL();
    void resizeGL(int w, int h);
    virtual void initializeGL();
    virtual void mousePressEvent(QMouseEvent *event);
    virtual void mouseReleaseEvent(QMouseEvent *event);
    virtual void mouseMoveEvent(QMouseEvent *event);
    virtual void wheelEvent(QWheelEvent * wheelEvent);
private:
    QTimer* m_timer = 0;
    bool m_doRender = false;
    glm::mat4 mat_mvp, mat_mv, mat_m, mat_v, mat_p;
    float m_aspectRatio = 1;
    std::shared_ptr<VertexBufferObject> m_vertices;
    std::shared_ptr<VertexIndexObject> m_indices;
    std::shared_ptr<VertexBufferObject> m_data;

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
    void receiveMesh(const MeshPackage & package);
};
#endif
