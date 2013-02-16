#include "../../../include/dec.hpp"
#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include "../../../include/render.hpp"
#include "../include/qtmainwindow.h"
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QFileDialog>
#include <QDir>
#include <QKeyEvent>
#include <iostream>
#include <random>
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
//#include <Eigen/ArpackSupport>
//#include <Eigen/CholmodSupport>
//#include "glutil.h"

MainWindow::MainWindow(QWidget * parent): QMainWindow(parent) {
    setMenuBar(new QMenuBar(this));
    QMenu *fileMenu = menuBar() -> addMenu(tr("&File"));
    QAction *openAct = new QAction(tr("&Open"), this);
    openAct->setShortcuts(QKeySequence::Open);
    openAct->setStatusTip(tr("Open"));
    connect(openAct, SIGNAL(triggered()), this, SLOT(openFile()));

    // Quit action
    QAction *quitAct = new QAction(tr("&Quit"), this);
    quitAct->setShortcuts(QKeySequence::Quit);
    quitAct->setStatusTip(tr("Quit"));
    connect(quitAct, SIGNAL(triggered()), this, SLOT(close()));

    fileMenu->addAction(openAct);
    fileMenu->addAction(quitAct);

    QGLFormat glFormat;
    std::cout << glFormat.openGLVersionFlags() << std::endl;
    std::cout << (QGLFormat::OpenGL_Version_3_0 <= glFormat.openGLVersionFlags()) << std::endl;
    if(QGLFormat::OpenGL_Version_3_3 & glFormat.openGLVersionFlags())
    {
        glFormat.setVersion( 3,3 );
    }
    else
    {
        glFormat.setVersion( 2, 1 );
    }
    std::cout << "GL Version: " << glFormat.majorVersion() << " " << glFormat.minorVersion() << std::endl;
    glFormat.setProfile( QGLFormat::CompatibilityProfile );
    glFormat.setSampleBuffers( true );
    m_glwidget = new GLWidget(glFormat,this);
    connect(this,SIGNAL(meshLoaded(const MeshPackage *))
            , m_glwidget,SLOT(recieveMesh(const MeshPackage *)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(recieveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);
}
#include <random>
//Eigen::MatrixXf m;
void MainWindow::openFile(const QString & filename) {
    MeshPackage *  package = new MeshPackage;
    m_mesh.reset(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    std::cout << "Read a mesh with verts:faces: " << m_mesh->vertices().size()<< ":" << m_mesh->numSimplices() << std::endl;
    m_dec.reset(new decltype(m_dec)::element_type(*m_mesh));
    typedef typename TriangleMeshf::Vector Vector;

    std::vector<Eigen::Vector3f> verts = m_mesh->vertices();
    mtao::normalizeInPlace(verts);//Normalize!!
    package = new MeshPackage{
        verts,
            mtao::template simplicesToRenderable<2>(*m_mesh),
            mtao::template simplicesToRenderable<1>(*m_mesh),
            verts,
            mtao::template simplicesToRenderable<2>(*m_mesh),
            mtao::template simplicesToRenderable<1>(*m_mesh)
    };







    m_glwidget->recieveMesh(package);
}






void MainWindow::openFile() {
    QFileDialog::Options options(QFileDialog::HideNameFilterDetails);
    QString filename = QFileDialog::getOpenFileName(
            this,
            tr("Choose file or directory"),
            QDir::homePath(),
            tr("OBJ Files (*.obj)")
            );

    if (filename.isNull()) {
        return;
    }
    openFile(filename);
}

void MainWindow::keyPressEvent(QKeyEvent *event) {
    switch(event->key()) {
        default: QMainWindow::keyPressEvent(event);return;
    }
    event->accept();
}
