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
    connect(this,SIGNAL(meshLoaded(std::shared_ptr<const MeshPackage>)), m_glwidget,SLOT(recieveMesh(std::shared_ptr<const MeshPackage>)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(recieveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);
}



void MainWindow::openFile(const QString & filename) {
    auto package = std::make_shared<MeshPackage>();
    m_mesh.reset(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    std::cout << "Read a mesh with verts:faces: " << m_mesh->vertices().size()<< ":" << m_mesh->numSimplices() << std::endl;
    m_dec.reset(new decltype(m_dec)::element_type(*m_mesh));
    typedef typename TriangleMeshf::Vector Vector;

    package->vertices = m_mesh->vertices();//copy here so we can normalize coordinates
    package->indices = mtao::template simplicesToRenderable<2>(*m_mesh);
    package->facevertices.resize(package->indices.size());
    package->faceindices.resize(package->indices.size());
    package->edgeindices = mtao::simplicesToRenderable<1>(*m_mesh);
    package->edgevertices.resize(package->edgeindices.size());



    int i;


    auto& verts = package->vertices;
    auto& packed_indices = package->indices;
    auto& faceverts = package->facevertices;
    auto& faceindices = package->faceindices;
    auto& edgeindices = package->edgeindices;
    auto& edgeverts = package->edgevertices;

    mtao::normalizeInPlace(verts);//Normalize!!

    std::transform(packed_indices.cbegin(), packed_indices.cend(), faceverts.begin(),
                   [&verts](const unsigned int ind)->typename std::remove_reference<decltype(faceverts)>::type::value_type
    {
        return verts[ind];
    });

    i=0;
    for(unsigned int & ind: faceindices)
    {
        ind = i++;
    }

    std::transform(edgeindices.cbegin(), edgeindices.cend(), edgeverts.begin(),
                   [&verts](const unsigned int ind)->typename std::remove_reference<decltype(edgeverts)>::type::value_type
    {
        return verts[ind];
    });

    i=0;
    for(unsigned int & ind: edgeindices)
    {
        ind = i++;
    }

    m_glwidget->recieveMesh(package);
    /*
    m_2form = m_dec->template genForm<PRIMAL_FORM,2>();
    m_2form.expr = decltype(m_2form.expr)::Zero(m_2form.expr.rows());
    m_1form = m_dec->template genForm<DUAL_FORM,1>();
    m_1form.expr = decltype(m_1form.expr)::Zero(m_1form.expr.rows());
    */
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

    switch(event->key()){
    default:
        QMainWindow::keyPressEvent(event);
        return;
    }
    event->accept();
}
