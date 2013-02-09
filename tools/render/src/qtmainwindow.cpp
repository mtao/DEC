#include "../include/qtmainwindow.h"
#include "../../../include/dec.hpp"
#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include "../../../include/render.hpp"
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QFileDialog>
#include <QDir>
#include <iostream>
#include <Eigen/SparseCholesky>
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
    m_glwidget = new GLWidget(this);
    connect(this,SIGNAL(meshLoaded(const MeshPackage &))
            , m_glwidget,SLOT(recieveMesh(const MeshPackage &)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(recieveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);
}
#include <random>
void MainWindow::openFile(const QString & filename) {
    MeshPackage package;
    std::unique_ptr<TriangleMeshf> mesh(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    DEC<TriangleMeshf> dec(*mesh);
    std::default_random_engine generator;
    std::uniform_real_distribution<float> distrib(-1,1);
    auto&& form2 = dec.template genForm<PRIMAL_FORM,2>();
    for(int i=0; i < form2.expr.rows(); ++i)//random data
        form2.expr(i) = distrib(generator);

    auto dhdh = dec.d(dec.h(dec.d(dec.template h<2>()))).expr.eval();
    typename Eigen::SimplicialLDLT<decltype(dhdh)> chol;
    chol.compute(dhdh);
    Eigen::VectorXf ret = chol.solve(form2.expr);//solve poisson problem
    form2.expr = ret;
    auto&& form1 = dec.template genForm<PRIMAL_FORM,1>();
    form1 = dec.h(dec.d(dec.h(form2)));//apply codifferential operator

    typedef typename TriangleMeshf::Vector Vector;
    const static int vecsize = sizeof(Vector)/sizeof(float);
    std::vector<Vector> verts = mesh->vertices();//copy here so we can normalize coordinates
    auto&& packed_indices = mtao::template simplicesToRenderable<2>(*mesh);


    mtao::normalizeInPlace(verts);//Normalize!!

    m_glwidget->makeCurrent();//activate glwidget opengl context for creating buffer objects
    package.vertices.reset(new VertexBufferObject((void*)verts.data(),
                                                  vecsize * verts.size(), GL_STATIC_DRAW));
    package.indices.reset(new VertexIndexObject((void*)packed_indices.data(), packed_indices.size(), GL_STATIC_DRAW, GL_TRIANGLES));


    int i;

    decltype(verts) faceverts(packed_indices.size());
    std::vector<unsigned int> faceindices(packed_indices.size());
    std::transform(packed_indices.cbegin(), packed_indices.cend(), faceverts.begin(),
                   [&verts](const unsigned int ind)->decltype(faceverts)::value_type
    {
        return verts[ind];
    });

    i=0;
    for(unsigned int & ind: faceindices)
    {
        ind = i++;
    }
    package.facevertices.reset(new VertexBufferObject((void*)faceverts.data(),
                                                      vecsize * faceverts.size(), GL_STATIC_DRAW));
    package.faceindices.reset(new VertexIndexObject((void*)faceindices.data(), faceindices.size(), GL_STATIC_DRAW, GL_TRIANGLES));
    auto&& packed_edge_indices = mtao::simplicesToRenderable<1>(*mesh);
    decltype(verts) edgeverts(packed_edge_indices.size());
    std::transform(packed_edge_indices.cbegin(), packed_edge_indices.cend(), edgeverts.begin(),
                   [&verts](const unsigned int ind)->decltype(edgeverts)::value_type
    {
        return verts[ind];
    });
    std::vector<unsigned int> edgeindices(edgeverts.size());

    i=0;
    for(unsigned int & ind: edgeindices)
    {
        ind = i++;
    }
    package.edgevertices.reset(new VertexBufferObject((void*)edgeverts.data(),
                                                      vecsize * edgeverts.size(), GL_STATIC_DRAW));
    package.edgeindices.reset(new VertexIndexObject((void*)edgeindices.data(), edgeindices.size(), GL_STATIC_DRAW, GL_LINES));

    emit meshLoaded(package);

    //dhdh.ldlt().solveInPlace(form2.expr);
    std::vector<std::array<float,3> > form2_ = mtao::formToRenderable(form2);

    FormPackage fpackage  = {"Test2",RT_FACE, std::make_shared<VertexBufferObject>(
                             (void*)form2_.data(),
                             3*form2_.size(),
                             GL_STATIC_DRAW,
                             1
                             )};
    emit formLoaded(fpackage);
    std::vector<std::array<float,2> > form1_ = mtao::formToRenderable(form1);
    fpackage  = {"Test1",RT_EDGE, std::make_shared<VertexBufferObject>(
                 form1_.data(),
                 2*form1_.size(),
                 GL_STATIC_DRAW,
                 1
                 )};
    emit formLoaded(fpackage);
    auto&& form0 = dec.template genForm<PRIMAL_FORM,0>();
    for(int i=0; i < form0.expr.rows(); ++i)
    {
        form0.expr(i) = -1;
    }
    std::vector<std::array<float,1> > form0_ = mtao::formToRenderable(form0);
    fpackage  = {"Test0",RT_VERT, std::make_shared<VertexBufferObject>(
                 form0_.data(),
                 1*form0_.size(),
                 GL_STATIC_DRAW,
                 1
                 )};
    emit formLoaded(fpackage);
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
    m_glwidget->keyPressEvent(event);
    QMainWindow::keyPressEvent(event);
}
