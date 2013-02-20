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

    //    mtao::normalizeInPlace(verts);//Normalize!!
    auto bbox = mtao::getBoundingBox(verts);
    mtao::normalizeToBBoxInPlace(verts,bbox);

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





    /*
    auto& dual_edgeindices = package->dual_edgeindices;
    auto& dual_edgeverts = package->dual_edgevertices;
    */
    auto& dual_vertices = package->dual_vertices;
    auto& dual_indices = package->dual_indices;
    auto& dual_edgeindices = package->dual_edgeindices;
    auto& dual_edgeverts = package->dual_edgevertices;
    auto& dual_faceindices = package->dual_faceindices;
    auto& dual_faceverts = package->dual_facevertices;//blarg! can't autosize because this might turn out to be stupidly sophisticated

    dual_vertices.resize(
                m_mesh->template numSimplices<2>() +
                m_mesh->template numSimplices<1>() +
                m_mesh->template numSimplices<0>()
                );
    int offset1 = m_mesh->template numSimplices<0>();
    int offset2 = offset1 + m_mesh->template numSimplices<1>();
    dual_edgeindices.resize(2*m_mesh->template numSimplices<1>());
    dual_edgeverts.resize(dual_edgeindices.size());


    typedef typename decltype(m_dec)::element_type::SparseMatrixColMajor SparseMatrix;
    const SparseMatrix & d0 = m_dec->template d<0>().expr;
    const SparseMatrix & d1 = m_dec->template d<1>().expr;
    /*
    const SparseMatrix & b2 = m_mesh->template b<2>();
    const SparseMatrix & b1 = m_mesh->template b<1>();
    */






    for(auto&& s: m_mesh->template simplices<0>()) {
        dual_vertices[s.Index()] = s.Center();
    }

    for(auto&& s: m_mesh->template simplices<1>()) {
        dual_vertices[s.Index()+offset1] = s.Center();
    }
    for(auto&& s: m_mesh->template simplices<2>()) {
        dual_vertices[s.Index()+offset2] = s.Center();
    }

    for(int i=0; i < m_mesh->template simplices<0>().size(); ++i) {
        auto&& s0 = m_mesh->template simplexByIndex<0>(i);

        for(SparseMatrix::InnerIterator it1(d0, s0.Index()); it1; ++it1) {
            auto&& s1 = m_mesh->template simplexByIndex<1>(it1.row());

            for(SparseMatrix::InnerIterator it2(d1, s1.Index()); it2; ++it2) {
                auto&& s2 = m_mesh->template simplexByIndex<2>(it2.row());
                dual_indices.push_back(s0.Index());
                dual_indices.push_back(s1.Index()+offset1);
                dual_indices.push_back(s2.Index()+offset2);
            }
        }
    }



    for(int i=0; i < m_mesh->template simplices<1>().size(); ++i) {
        auto&& s = m_mesh->template simplexByIndex<1>(i);
        SparseMatrix::InnerIterator it(d1, s.Index());
        auto&& s1 = m_mesh->template simplexByIndex<2>(it.row());
        dual_edgeindices[2*i] = 2*i;
        dual_edgeverts[2*i] = s1.Center();
        ++it;
        auto&& s2 = m_mesh->template simplexByIndex<2>(it.row());
        dual_edgeindices[2*i+1] = 2*i+1;
        dual_edgeverts[2*i+1] = s2.Center();

    }







    dual_faceverts.clear();
    dual_faceverts.reserve(3*m_mesh->template numSimplices<2>());
    dual_faceindices.clear();
    dual_faceverts.reserve(3*m_mesh->template numSimplices<2>());

    m_dual_vertex_form_indices.clear();
    m_dual_vertex_form_indices.reserve(m_mesh->template numSimplices<0>());

    int count=0;
    for(int i=0; i < m_mesh->template simplices<0>().size(); ++i) {
        auto&& s0 = m_mesh->template simplexByIndex<0>(i);

        for(SparseMatrix::InnerIterator it1(d0, s0.Index()); it1; ++it1) {
            auto&& s1 = m_mesh->template simplexByIndex<1>(it1.row());

            for(SparseMatrix::InnerIterator it2(d1, s1.Index()); it2; ++it2) {
                auto&& s2 = m_mesh->template simplexByIndex<2>(it2.row());

                dual_faceverts.push_back(s0.Center());
                dual_faceverts.push_back(s1.Center());
                dual_faceverts.push_back(s2.Center());
                ++count;
            }
        }
        m_dual_vertex_form_indices.push_back(count);
    }


    dual_faceindices.clear();
    dual_faceindices.resize(dual_faceverts.size());
    for(int i=0; i < dual_faceindices.size(); ++i) {
        dual_faceindices[i] = i;
    }





    mtao::normalizeToBBoxInPlace(dual_faceverts,bbox);
    mtao::normalizeToBBoxInPlace(dual_edgeverts,bbox);
    mtao::normalizeToBBoxInPlace(dual_vertices,bbox);





























    emit meshLoaded(package);
    auto m_2form = m_dec->template genForm<PRIMAL_FORM,2>();
    m_2form.expr = decltype(m_2form.expr)::Zero(m_2form.expr.rows());
    auto m_1form = m_dec->template genForm<DUAL_FORM,1>();
    m_1form.expr = decltype(m_1form.expr)::Zero(m_1form.expr.rows());
    m_1form.expr = decltype(m_1form.expr)::Ones(m_1form.expr.rows());
    auto m_0form = m_dec->template genForm<DUAL_FORM,0>();
    m_0form.expr = decltype(m_0form.expr)::Zero(m_0form.expr.rows());
    emit formLoaded(makeFormPackage("test0", m_0form));
    emit formLoaded(makeFormPackage("test1", m_1form));
    emit formLoaded(makeFormPackage("test2", m_2form));
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
        m_glwidget->keyPressEvent(event);
        QMainWindow::keyPressEvent(event);
        return;
    }
    event->accept();
}
