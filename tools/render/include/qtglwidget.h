#ifndef GLWIDGET_H
#define GLWIDGET_H
#include "glutil.h"
#include <QGLWidget>
#include <glm/glm.hpp>
#include <memory>
#include <QTimer>
#include <QTime>
#include <QVector3D>

#include "packages.h"
bool operator<(const FormPackage & a, const FormPackage & b);

class MainWindow;
class GLWidget: public QGLWidget
{
    Q_OBJECT
public:
    friend class MainWindow;
    GLWidget(const QGLFormat&  format, QWidget * parent=0);
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
    //void render(RenderType type);
    void renderForm(const FormPackage & form);
    QTimer* m_timer = 0;
    int m_renderType = RT_NONE;
    bool m_doRender = false;
    glm::mat4 mat_mvp, mat_mv, mat_m, mat_v, mat_p;
    float m_aspectRatio = 1;
    std::unique_ptr<ShaderProgram> m_shader;
    std::unique_ptr<ShaderProgram> m_vertshader;
    std::unique_ptr<ShaderProgram> m_faceshader;
    std::unique_ptr<ShaderProgram> m_edgeshader;
    MeshBuffers m_meshbuffers;
    std::shared_ptr<const MeshPackage> m_meshpackage;
    std::map<QString, FormPackage> m_formpackages;
    std::set<QString> m_active_forms;

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
    void receiveMesh(std::shared_ptr<const MeshPackage> package);
    void receiveForm(const FormPackage & package);
    void enableForm(const QString & formname);
    void disableForm(const QString & name);
    void clearForms();
};
#endif
