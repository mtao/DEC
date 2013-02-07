#include "../include/glwidget.h"
#include <glm/gtc/matrix_transform.hpp>
#include <QMouseEvent>
#include <qmath.h>

GLWidget::GLWidget(QWidget * parent): QGLWidget(parent), m_timer(new QTimer(this)) {

    connect( m_timer, SIGNAL( timeout() ), SLOT( update() ) );

    m_timer->start( 16 );
}
void GLWidget::paintGL() {
    glClearColor(0.0,0.0,0.0,0.0);
    glClear( GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT );
    const int delta = m_time.elapsed() - m_lastTime;
    m_rotation += m_angularMomentum * (delta / 1000.0);
    m_lastTime += delta;

    // Setup the modelview matrix
    mat_m = glm::translate(glm::mat4(1.0f),glm::vec3(0.f,0.f,-m_distance));
    mat_m = glm::rotate(mat_m, (float)m_rotation.x(), glm::vec3(1.f,0.f,0.f));
    mat_m = glm::rotate(mat_m, (float)m_rotation.y(), glm::vec3(0.f,1.f,0.f));
    mat_m = glm::rotate(mat_m, (float)m_rotation.z(), glm::vec3(0.f,0.f,1.f));
}

void GLWidget::resizeGL(int w, int h) {
    // Set the viewport to window dimensions
    glViewport( 0, 0, w, h );

    // Setup an orthogonal projection matrix
    float viewingAngle = 40.0;
    float nearPlane = 0.0001;
    float farPlane = 100.0;
    h = std::max( h, 1 );
    m_aspectRatio = float( w ) / float( h );
    mat_p = glm::perspective( viewingAngle, m_aspectRatio, nearPlane, farPlane );
    mat_mvp = mat_p * mat_mv;
}

void GLWidget::initializeGL() {
    glm::vec3 eyePosition(0.0f,0.0f,1.0f);
    glm::vec3 targetPosition(0.0f,0.0f,0.0f);
    glm::vec3 upDirection(0.0f,1.0f,0.0f);
    mat_v = glm::lookAt(eyePosition,targetPosition,upDirection);
}
void GLWidget::mouseMoveEvent(QMouseEvent *event)
{
    if (event->buttons() & Qt::LeftButton) {
        const QPointF delta = event->posF() - lastPos;
        lastPos=event->posF();
        const QVector3D angularImpulse = QVector3D(delta.y(), delta.x(), 0) * 0.1;

        m_rotation += angularImpulse;
        m_accumulatedMomentum += angularImpulse;

        event->accept();
        update();
    }
}

void GLWidget::mousePressEvent(QMouseEvent *event)
{

    m_mouseEventTime = m_time.elapsed();
    m_angularMomentum = m_accumulatedMomentum = QVector3D();
    event->accept();
    lastPos=event->posF();
}

void GLWidget::mouseReleaseEvent(QMouseEvent *event)
{

    const int delta = m_time.elapsed() - m_mouseEventTime;
    m_angularMomentum = m_accumulatedMomentum * (1000.0 / qMax(1, delta));
    event->accept();
    update();
}

void GLWidget::wheelEvent(QWheelEvent *event)
{

    m_distance *= qPow(1.2, -event->delta() / 120);
    event->accept();
    update();
}

void GLWidget::receiveMesh(const MeshPackage &package) {
    m_vertices = package.vertices;
    m_indices = package.indices;
}
