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
#include <Eigen/ArpackSupport>
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
    connect(this,SIGNAL(meshLoaded(const MeshPackage &))
            , m_glwidget,SLOT(recieveMesh(const MeshPackage &)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(recieveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);
}
#include <random>
//Eigen::MatrixXf m;
void MainWindow::openFile(const QString & filename) {
    MeshPackage package;
    m_mesh.reset(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    std::cout << "Read a mesh with verts:faces: " << m_mesh->vertices().size()<< ":" << m_mesh->numSimplices() << std::endl;
    //m_dec.reset(new DEC<TriangleMeshf>(*m_mesh));
    m_dec.reset(new decltype(m_dec)::element_type(*m_mesh));
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
    m_1form = m_dec->template genForm<DUAL_FORM,1>();
    m_1form.expr = decltype(m_1form.expr)::Zero(m_1form.expr.rows());
    /*
    for(auto&& s: m_mesh->template simplices<1>()) {
        std::cout << s.Volume() << " " << s.DualVolume() << std::endl;
    }
    */

    /*
    Eigen::SparseMatrix<float> dhdh = m_dec->d(m_dec->h(m_dec->d(m_dec->template h<2>()))).expr.eval();
    Eigen::ArpackGeneralizedSelfAdjointEigenSolver<
            Eigen::SparseMatrix<float>,
            typename Eigen::SimplicialLDLT<decltype(dhdh)>,
            true
            > eigensolver(dhdh, 100);
    std::cout << eigensolver.eigenvalues().transpose() << std::endl;
    m = eigensolver.eigenvectors();

    */
}





#define CHECK_INTERIORS

bool dumb = true;
void MainWindow::randomData() {
    if(!m_dec || !m_mesh) {return;}
    auto&& form2 = m_dec->template genForm<PRIMAL_FORM,2>();
    auto&& form1 = m_dec->template genForm<PRIMAL_FORM,1>();
    auto&& form0 = m_dec->template genForm<PRIMAL_FORM,0>();
#ifdef CHECK_INTERIORS
    auto&& interior2 = m_mesh->template interior<2>();
    for(int i=0; i < interior2.rows(); ++i) {

        form2.expr(i) = interior2.diagonal()(i);
        if(dumb)
            form2.expr(i) = 1-form2.expr(i);
    }
    auto&& interior1 = m_mesh->template interior<1>();
    for(int i=0; i < interior1.rows(); ++i) {
        form1.expr(i) = interior1.diagonal()(i);
        if(dumb)
            form1.expr(i) = 1-form1.expr(i);
    }
    auto&& interior0 = m_mesh->template interior<0>();
    for(int i=0; i < interior0.rows(); ++i) {
        form0.expr(i) = interior0.diagonal()(i);
        if(dumb)
            form0.expr(i) = 1-form0.expr(i);
    }
    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test0",form0));
    return;
#endif
    //form2.expr = decltype(form2.expr)::Random(form2.expr.rows());

    Eigen::SparseMatrix<float> dhdh = m_dec->h(m_dec->d(m_dec->h(m_dec->template d<1>()))).expr.eval();
    /*
    if(dumb) {
        form2 = m_2form;
        dumb = false;
    std::cout << dhdh << std::endl;
    std::cout << m_dec->template h<2>().expr << std::endl;
    std::cout << m_dec->template d<0,DUAL_FORM>().expr << std::endl;
    std::cout << m_dec->template h<1,DUAL_FORM>().expr << std::endl;
    std::cout << m_dec->template d<1,PRIMAL_FORM>().expr << std::endl;
    for(auto&& s: m_mesh->template simplices<1>()) {
        std::cout << s.Volume() << " " << s.DualVolume() << std::endl;
    }
    }
    */

        /*
    form2.expr = dhdh*form2.expr;
    std::cout << form2.expr.norm() << std::endl;
    form2.expr.normalize();

    */



    //typename Eigen::SimplicialLDLT<decltype(dhdh)> chol;
    typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::SimplicialLDLT<decltype(dhdh)> > chol;
    chol.setMaxIterations(dhdh.rows()/5);
    chol.setTolerance(0.1);
    //typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::CholmodBaseSupernodalLLT<decltype(dhdh)> > chol;
    chol.compute(dhdh);
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at dec->mposition" << std::endl;
    }
    Eigen::VectorXf ret = chol.solve(m_1form.expr);//solve poisson problem
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at solving" << std::endl;
    }
    m_1form.expr = ret;
    //std::cout << "Norm error: " << (dhdh* ret - m_1form.expr).norm() << std::endl;
    form1.expr = ret / ret.lpNorm<Eigen::Infinity>();
    /*
    static int WHICH = 0;
    form2.expr = m.col(WHICH);
    std::cout << form2.expr.transpose() << std::endl;
    WHICH = (WHICH+1)%100;
    */
    /*
    form2 = m_dec->d(m_dec->h(m_dec->d(m_dec->h(m_2form))));
    form2.expr = form2.expr / form2.expr.lpNorm<Eigen::Infinity>();
    */
    //form2 = m_dec->d(m_dec->h(form1));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));
    return;
    m_2form = form2;

    //form1 = m_dec->h(m_dec->d(m_dec->h(form2)));//apply codifferential operator
    //m_2form = form2;
    //form2 = m_2form;

    for(int i=0; i < form1.expr.rows(); ++i) {
        form1.expr(i) = m_mesh->template simplexByIndex<1>(i).DualVolume();
    }
    form1.expr = form1.expr / form1.expr.lpNorm<Eigen::Infinity>();

    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));

    for(int i=0; i < form0.expr.rows(); ++i)
    {
        form0.expr(i) = -1;
    }
    m_glwidget->recieveForm(mtao::makeFormPackage("Test0",form0));
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
        std::uniform_int_distribution<int> rand;
    switch(event->key()) {
    case Qt::Key_R:
        randomData();
        break;
    case Qt::Key_T:
        m_1form.expr = decltype(m_1form.expr)::Zero(m_1form.expr.rows());
        rand = std::uniform_int_distribution<int>(0,m_1form.expr.rows()-1);
        //for(int i=0; i < rand(generator); ++i)
        m_1form.expr(rand(generator)) = 10;
    //m_2form.expr = decltype(m_2form.expr)::Random(m_2form.expr.rows());
#ifdef CHECK_INTERIORS
        dumb ^= true;
        randomData();
        return;
#else
        dumb = true;
#endif


    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",m_1form));
        break;

    default:
        m_glwidget->keyPressEvent(event);
        QMainWindow::keyPressEvent(event);
        return;
    }
    event->accept();
}
