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
#include <QKeyEvent>
#include <iostream>
#include <random>
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
    connect(this,SIGNAL(meshLoaded(const MeshPackage &))
            , m_glwidget,SLOT(recieveMesh(const MeshPackage &)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(recieveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);
}
#include <random>
void MainWindow::openFile(const QString & filename) {
    MeshPackage package;
    m_mesh.reset(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    std::cout << "Read a mesh with verts:faces: " << m_mesh->vertices().size()<< ":" << m_mesh->numSimplices() << std::endl;
    m_dec.reset(new DEC<TriangleMeshf>(*m_mesh));
    typedef typename TriangleMeshf::Vector Vector;
    const static int vecsize = sizeof(Vector)/sizeof(float);
    std::vector<Vector> verts = m_mesh->vertices();//copy here so we can normalize coordinates
    auto&& packed_indices = mtao::template simplicesToRenderable<2>(*m_mesh);


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
    auto&& packed_edge_indices = mtao::simplicesToRenderable<1>(*m_mesh);
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

    m_glwidget->recieveMesh(package);
    m_2form = m_dec->template genForm<PRIMAL_FORM,2>();
    m_2form.expr = decltype(m_2form.expr)::Zero(m_2form.expr.rows());
}
void MainWindow::randomData() {
    if(!m_dec || !m_mesh) {return;}
    auto&& form2 = m_dec->template genForm<PRIMAL_FORM,2>();
    form2.expr = decltype(form2.expr)::Random(form2.expr.rows());

    Eigen::SparseMatrix<float> dhdh = m_dec->d(m_dec->h(m_dec->d(m_dec->template h<2>()))).expr.eval();





    typename Eigen::SimplicialLDLT<decltype(dhdh)> chol;
    chol.compute(dhdh);
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at dec->mposition" << std::endl;
    }
    //Eigen::VectorXf ret = chol.solve(form2.expr);//solve poisson problem
    Eigen::VectorXf ret = chol.solve(m_2form.expr);//solve poisson problem
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at solving" << std::endl;
    }
    form2.expr = ret;
    auto&& form1 = m_dec->template genForm<PRIMAL_FORM,1>();
    form1 = m_dec->h(m_dec->d(m_dec->h(form2)));//apply codifferential operator
    form1.expr /= form1.expr.template lpNorm<Eigen::Infinity>();
    //form2.expr/= (form2.expr.template lpNorm<Eigen::Infinity>()) /2;
    m_2form = form2;
    //m_2form.expr = dhdh * m_2form.expr;

    //form2.expr/= (form2.expr.template lpNorm<Eigen::Infinity>());
    //form2 = m_2form;


    //dhdh.ldlt().solveInPlace(form2.expr);
    std::vector<std::array<float,3> > form2_ = mtao::formToRenderable(form2);

    FormPackage fpackage  = {"Test2",RT_FACE, std::make_shared<VertexBufferObject>(
                             (void*)form2_.data(),
                             3*form2_.size(),
                             GL_STATIC_DRAW,
                             1
                             )};
    m_glwidget->recieveForm(fpackage);
    std::vector<std::array<float,2> > form1_ = mtao::formToRenderable(form1);
    fpackage  = {"Test1",RT_EDGE, std::make_shared<VertexBufferObject>(
                 form1_.data(),
                 2*form1_.size(),
                 GL_STATIC_DRAW,
                 1
                 )};
    m_glwidget->recieveForm(fpackage);
    auto&& form0 = m_dec->template genForm<PRIMAL_FORM,0>();
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
    m_glwidget->recieveForm(fpackage);
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
    static std::default_random_engine generator;
    FormPackage fpackage;
        std::vector<std::array<float,3> > m_2form_;
        std::uniform_int_distribution<int> rand;
    switch(event->key()) {
    case Qt::Key_R:
        randomData();
        break;
    case Qt::Key_T:
        m_2form.expr = decltype(m_2form.expr)::Zero(m_2form.expr.rows());
        rand = std::uniform_int_distribution<int>(0,m_2form.expr.rows());
        for(int i=0; i < rand(generator); ++i)
        m_2form.expr(rand(generator)) = 1;
        m_2form_ = mtao::formToRenderable(m_2form);

        fpackage  = {"Test2",RT_FACE, std::make_shared<VertexBufferObject>(
                                 (void*)m_2form_.data(),
                                 3*m_2form_.size(),
                                 GL_STATIC_DRAW,
                                 1
                                 )};
        m_glwidget->recieveForm(fpackage);
        break;

    default:
        m_glwidget->keyPressEvent(event);
        QMainWindow::keyPressEvent(event);
        return;
    }
    event->accept();
}