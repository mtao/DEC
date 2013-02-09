#include "../include/qtglwidget.h"
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <QMouseEvent>
#include <QDebug>
#include <QResource>
#include <memory>
#include <utility>
#include <qmath.h>
#include <QFileInfo>
#include <QDir>

GLWidget::GLWidget(const QGLFormat & format, QWidget * parent): QGLWidget(format, parent), m_timer(new QTimer(this)) {

    connect( m_timer, SIGNAL( timeout() ), SLOT( update() ) );

    m_timer->start( 16 );
    QResource::registerResource("../qrc_shaders.cxx");
    QDir d(":/");

}

void GLWidget::resizeGL(int w, int h) {
    // Set the viewport to window dimensions
    glViewport( 0, 0, w, h );

    // Setup an orthogonal projection matrix
    float viewingAngle = 40.0;
    float nearPlane = 0.1;
    float farPlane = 10000.0;
    h = std::max( h, 1 );
    m_aspectRatio = float( w ) / float( h );
    mat_p = glm::perspective( viewingAngle, m_aspectRatio, nearPlane, farPlane );
    //mat_mvp = mat_p * mat_mv;
}

void GLWidget::initializeGL() {
    glewInit();
    glm::vec3 eyePosition(0.0f,0.0f,1.0f);
    glm::vec3 targetPosition(0.0f,0.0f,0.0f);
    glm::vec3 upDirection(0.0f,1.0f,0.0f);
    mat_v = glm::lookAt(eyePosition,targetPosition,upDirection);


    m_shader.reset(new ShaderProgram());
    m_faceshader.reset(new ShaderProgram());
    m_vertshader.reset(new ShaderProgram());
    m_edgeshader.reset(new ShaderProgram());
    initShader(*m_shader, "");
    initShader(*m_faceshader, "face");
    initShader(*m_vertshader, "vert");
    initShader(*m_edgeshader, "edge");

    glEnable(GL_DEPTH_TEST);
    glDisable(GL_CULL_FACE);
}
void GLWidget::initShader(ShaderProgram & program, const QString & geotype)
{
    QString vertexShaderPath(":/shader.v.glsl");
    QString fragmentShaderPath(":/shader.f.glsl");
    QString geometryShaderPath(tr(":/") + geotype + tr("shader.g.glsl"));
    if(this->format().majorVersion() < 2 || (geotype != tr("") && geotype != tr("face"))) {
        vertexShaderPath = tr(":/shader.130.v.glsl");
        fragmentShaderPath = tr(":/shader.130.f.glsl");
        geometryShaderPath = tr("");//hopefully "" file won't exist
    }
    GLuint programId = program.programId;
    // First we load and compile the vertex shader...
    GLuint vertexId = compileShader(GL_VERTEX_SHADER, vertexShaderPath);
    // ...now the fragment shader...
    GLuint fragmentId = compileShader(GL_FRAGMENT_SHADER, fragmentShaderPath );


    glUseProgram(programId);
    if(vertexId && fragmentId)
    {
        glAttachShader(programId, vertexId);
        glAttachShader(programId, fragmentId);
    }
    else
    {
        glDeleteProgram(programId);
        return;
    }

    QFileInfo geoFileInfo(geometryShaderPath);
    GLuint geometryId = 0;
    if(geoFileInfo.exists())
    {
        geometryId = compileShader(GL_GEOMETRY_SHADER, geometryShaderPath );
        if(geometryId)
            glAttachShader(programId, geometryId);
    }
    glLinkProgram(programId);
    GLint link_ok = GL_FALSE;
    glGetProgramiv(programId, GL_LINK_STATUS, &link_ok);

    if(link_ok == GL_FALSE)
    {
        GLint log_length = 0;
        glGetProgramiv(programId, GL_INFO_LOG_LENGTH, &log_length);
        char * log = new char[log_length];
        glGetProgramInfoLog(programId, log_length, NULL, log);
        qWarning() << log;
        delete[] log;
        glDeleteProgram(programId);
        return;
    }

    glUseProgram(0);





}


GLuint GLWidget::compileShader(GLenum shaderType, const QString & fileName)
{

    qWarning() << fileName ;
    QFile file(fileName);
    if(!file.open(QIODevice::ReadOnly | QIODevice::Text))
    {
        qWarning() << "File not found: " << fileName;
        return false;
    }
    QString contents = file.readAll();
    GLuint shaderId = glCreateShader(shaderType);
    char * buf = new char[contents.length()+1];
    memcpy(buf,contents.toStdString().c_str(),contents.length());
    buf[contents.length()] = '\0';
    glShaderSource(shaderId,1, const_cast<const char **>(&buf),NULL);
    delete[] buf;
    glCompileShader(shaderId);
    GLint compile_ok;
    glGetShaderiv(shaderId, GL_COMPILE_STATUS, &compile_ok);
    if(compile_ok == GL_FALSE)
    {
        GLint log_length = 0;
        glGetShaderiv(shaderId, GL_INFO_LOG_LENGTH, &log_length);
        char * log = new char[log_length];
        glGetShaderInfoLog(shaderId, log_length, NULL, log);
        qWarning() << log;
        delete[] log;
    }

    return shaderId;
}



void GLWidget::recieveMesh(const MeshPackage &package) {
    m_meshpackage = package;
    m_shader->addAttribute("vertex", package.vertices);
    m_vertshader->addAttribute("vertex", package.vertices);
    m_edgeshader->addAttribute("vertex", package.edgevertices);
    m_faceshader->addAttribute("vertex", package.facevertices);
}


void GLWidget::recieveForm(const FormPackage &package) {
    m_formpackages[package.title] = package;
    enableForm(package.title);//TODO: make an interface for enabling/disabling forms instead of rendering the last n-form we've seen
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
    mat_mv = mat_v * mat_m;
    mat_mvp = mat_p * mat_mv;


    if(!m_meshpackage.indices) return;
    /*

    ShaderProgram * shader;
    switch(m_renderType)
    {
    case RT_FACE: shader = m_faceshader.get(); break;
    case RT_EDGE: shader = m_edgeshader.get(); break;
    case RT_VERT: shader = m_vertshader.get(); break;
    }
    shader->bind();
    glUniformMatrix4fv(glGetUniformLocation(shader->programId, "MVP"),
                       1, GL_FALSE, glm::value_ptr(mat_mvp));

    switch(m_renderType)
    {
    case RT_FACE: m_meshpackage.faceindices->render(); break;
    case RT_EDGE: m_meshpackage.edgeindices->render(); break;
    case RT_VERT: glDrawArrays(GL_POINTS, 0, m_meshpackage.vertices->size); break;
    }
    shader->release();
    */
    for(RenderType t: {RT_FACE, RT_EDGE, RT_VERT}) {
        if(m_renderType & t) {
            render(t);
        }
    }
    if(m_renderType == 0) {
        render(RT_NONE);
    }

}

std::unique_ptr<ShaderProgram> & GLWidget::shaderSelector(RenderType type) {

    switch(type)
    {
    case RT_FACE: return m_faceshader;
    case RT_EDGE: return m_edgeshader;
    case RT_VERT: return m_vertshader;
    case RT_NONE: return m_shader;
    }
}

/*
void GLWidget::renderForm(const FormPackage & form) {

    ShaderProgram * shader = shaderSelector(form.type);
    shader->bind();
    glUniformMatrix4fv(glGetUniformLocation(shader->programId, "MVP"),
                       1, GL_FALSE, glm::value_ptr(mat_mvp));
    GLint attributeId = glGetAttribLocation(shader->programId, "data");
    form.data->bind(attributeId);

    switch(form.type)
    {
    case RT_FACE: m_meshpackage.faceindices->render(); break;
    case RT_EDGE: m_meshpackage.edgeindices->render(); break;
    case RT_VERT: glDrawArrays(GL_POINTS, 0, m_meshpackage.vertices->size); break;
    }
    shader->release();
}
*/
void GLWidget::render(RenderType type) {

    auto&& shader = shaderSelector(type);
    shader->bind();
    glUniformMatrix4fv(glGetUniformLocation(shader->programId, "MV"),
                       1, GL_FALSE, glm::value_ptr(mat_mv));
    glUniformMatrix4fv(glGetUniformLocation(shader->programId, "MVP"),
                       1, GL_FALSE, glm::value_ptr(mat_mvp));
    switch(type)
    {
    case RT_NONE: m_meshpackage.indices->render(); break;
    case RT_FACE: m_meshpackage.faceindices->render(); break;
    case RT_EDGE: m_meshpackage.edgeindices->render(); break;
    case RT_VERT: glDrawArrays(GL_POINTS, 0, m_meshpackage.vertices->size); break;
    }
    shader->release();
}





void GLWidget::keyPressEvent(QKeyEvent *event) {
    switch(event->key()) {
    case Qt::Key_F: m_renderType ^= RT_FACE; break;
    case Qt::Key_E: m_renderType ^= RT_EDGE; break;
    case Qt::Key_V: m_renderType ^= RT_VERT; break;
    default: return;
    }
    //If I haven't used the key I'll have already returned, so I should accept the input
    event->accept();
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

void GLWidget::enableForm(const QString &formname) {
    auto&& form = m_formpackages[formname];
    shaderSelector(form.type)->addAttribute("data", form.data);
    //m_renderType ^= form.type;
}
