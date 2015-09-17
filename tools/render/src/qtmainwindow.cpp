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
#include <QDockWidget>
#include <iostream>
#include <random>

MainWindow::MainWindow(QWidget * parent, FormBar * bar): QMainWindow(parent), m_formbar(bar) {
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
    connect(this,SIGNAL(meshLoaded(std::shared_ptr<const MeshPackage>)), m_glwidget,SLOT(receiveMesh(std::shared_ptr<const MeshPackage>)));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_glwidget,SLOT(receiveForm(const FormPackage &)));
    setCentralWidget(m_glwidget);



    QDockWidget * dock = new QDockWidget(tr("Form Chooser"),this);
    if(!m_formbar) {
        m_formbar = new FormBar(this);
    }
    dock->setWidget(m_formbar);
    addDockWidget(Qt::LeftDockWidgetArea, dock);
    connect(this,SIGNAL(loadingNewMesh())
            , m_glwidget,SLOT(unloadMesh()));
    connect(this,SIGNAL(formLoaded(const FormPackage &))
            , m_formbar,SLOT(receiveForm(const FormPackage &)));
    connect(this,SIGNAL(particlesLoaded(std::shared_ptr<VertexBufferObject>))
            , m_glwidget,SLOT(receiveParticles(std::shared_ptr<VertexBufferObject>)));
    connect(
            m_formbar, SIGNAL(enableForm(const QString &)),
            m_glwidget, SLOT(enableForm(const QString &)));
    connect(
            m_formbar, SIGNAL(disableForm(const QString &)),
            m_glwidget, SLOT(disableForm(const QString &)));
    connect(
            m_formbar, SIGNAL(clearForms(void)),
            m_glwidget, SLOT(clearForms(void)));

}



void MainWindow::openFile(const QString & filename) {
    emit loadingNewMesh();
    m_mesh.reset(readOBJtoSimplicialComplex<MeshType>(filename.toStdString()));
    std::cout << "Read a mesh with verts:faces: " << m_mesh->vertices().size()<< ":" << m_mesh->numSimplices() << std::endl;
    initializeMesh();
}

void MainWindow::initializeMesh() {
    auto package = std::make_shared<MeshPackage>();
    m_dec.reset(new decltype(m_dec)::element_type(*m_mesh));
    typedef typename MeshType::Vector Vector;

    package->vertices = m_mesh->vertices();//copy here so we can normalize coordinates
    package->indices = mtao_internal::template simplicesToRenderable<2>(*m_mesh);
    package->facevertices.resize(package->indices.size());
    package->faceindices.resize(package->indices.size());
    package->edgeindices = mtao_internal::simplicesToRenderable<1>(*m_mesh);
    package->edgevertices.resize(package->edgeindices.size());



    int i;


    auto& verts = package->vertices;
    auto& packed_indices = package->indices;
    auto& faceverts = package->facevertices;
    auto& faceindices = package->faceindices;
    auto& edgeindices = package->edgeindices;
    auto& edgeverts = package->edgevertices;

    //    mtao::normalizeInPlace(verts);//Normalize!!
    std::cout << "Normalizing vertices..." << std::endl;
    m_bbox = mtao::getBoundingBox(verts);
    mtao::normalizeToBBoxInPlace(verts,m_bbox);

    std::cout << "Writing primal verts" << std::endl;
    std::transform(packed_indices.cbegin(), packed_indices.cend(), faceverts.begin(),
                   [&verts](const unsigned int ind)->typename std::remove_reference<decltype(faceverts)>::type::value_type
    {
        return verts[ind];
    });

    i=0;
    std::cout << "Writing primal face indices" << std::endl;
    for(unsigned int & ind: faceindices)
    {
        ind = i++;
    }

    std::cout << "Writing primal edge vertices" << std::endl;
    std::transform(edgeindices.cbegin(), edgeindices.cend(), edgeverts.begin(),
                   [&verts](const unsigned int ind)->typename std::remove_reference<decltype(edgeverts)>::type::value_type
    {
        return verts[ind];
    });

    i=0;
    std::cout << "Writing primal edge indices" << std::endl;
    for(unsigned int & ind: edgeindices)
    {
        ind = i++;
    }





    std::cout << "Writing dual quantities" << std::endl;
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
    int offset2 = 0;//offset1 + m_mesh->template numSimplices<1>();
    int offset1 = m_mesh->template numSimplices<2>();
    int offset0 = offset1 + m_mesh->template numSimplices<1>();
    dual_edgeindices.resize(2*m_mesh->template numSimplices<1>());
    dual_edgeverts.resize(dual_edgeindices.size());


    typedef typename decltype(m_dec)::element_type::SparseMatrixColMajor SparseMatrix;
    const SparseMatrix & d0 = m_dec->template d<0>().expr;
    const SparseMatrix & d1 = m_dec->template d<1>().expr;





    package->num_dual_verts = m_mesh->template numSimplices<2>();

    std::cout << "Writing dual verts" << std::endl;
    for(auto&& s: m_mesh->template simplices<2>()) {
        dual_vertices[s.Index()+offset2] = s.Center();
    }

    for(auto&& s: m_mesh->template simplices<1>()) {
        dual_vertices[s.Index()+offset1] = s.Center();
    }

    for(auto&& s: m_mesh->template simplices<0>()) {
        dual_vertices[s.Index()+offset0] = s.Center();
    }

    std::cout << "Writing dual mesh indices" << std::endl;
    for(int i=0; i < m_mesh->template simplices<0>().size(); ++i) {
        auto&& s0 = m_mesh->template simplex<0>(i);

        for(SparseMatrix::InnerIterator it1(d0, s0.Index()); it1; ++it1) {
            auto&& s1 = m_mesh->template simplex<1>(it1.row());

            for(SparseMatrix::InnerIterator it2(d1, s1.Index()); it2; ++it2) {
                auto&& s2 = m_mesh->template simplex<2>(it2.row());
                dual_indices.push_back(s0.Index()+offset0);
                dual_indices.push_back(s1.Index()+offset1);
                dual_indices.push_back(s2.Index()+offset2);
            }
        }
    }



    std::cout << "Writing dual edge indices" << std::endl;
    for(int i=0; i < m_mesh->template simplices<1>().size(); ++i) {
        auto&& s = m_mesh->template simplex<1>(i);
        SparseMatrix::InnerIterator it(d1, s.Index());
        if(!it) continue;
        auto&& s1 = m_mesh->template simplex<2>(it.row());
        dual_edgeindices[2*i] = 2*i;
        dual_edgeverts[2*i] = s1.Center();
        ++it;
        if(!it) continue;
        auto&& s2 = m_mesh->template simplex<2>(it.row());
        dual_edgeindices[2*i+1] = 2*i+1;
        dual_edgeverts[2*i+1] = s2.Center();

    }






    std::cout << "Reserving" << std::endl;

    dual_faceverts.clear();
    dual_faceverts.reserve(3*m_mesh->template numSimplices<2>());
    dual_faceindices.clear();
    dual_faceindices.reserve(3*m_mesh->template numSimplices<2>());

    m_dual_vertex_form_indices.clear();
    m_dual_vertex_form_indices.reserve(m_mesh->template numSimplices<0>());

    int count=0;
    std::cout << "Writing dual face indices" << std::endl;
    for(int i=0; i < m_mesh->template simplices<0>().size(); ++i) {
        auto&& s0 = m_mesh->template simplex<0>(i);

        for(SparseMatrix::InnerIterator it1(d0, s0.Index()); it1; ++it1) {
            auto&& s1 = m_mesh->template simplex<1>(it1.row());

            SparseMatrix::InnerIterator it2(d1, s1.Index());
        if(!it2) continue;

            auto&& s2 = m_mesh->template simplex<2>(it2.row());
            dual_faceverts.push_back(s0.Center());
            if((s0[0] == s1[0]) ^ s1.isSameSign(s2)) {
            dual_faceverts.push_back(s2.Center());
            dual_faceverts.push_back(s1.Center());
            } else {
            dual_faceverts.push_back(s1.Center());
            dual_faceverts.push_back(s2.Center());
            }
            ++count;
            ++it2;
        if(!it2) continue;

            auto&& s2b = m_mesh->template simplex<2>(it2.row());
            dual_faceverts.push_back(s0.Center());
            if((s0[0] == s1[0]) ^ s1.isSameSign(s2)) {
            dual_faceverts.push_back(s1.Center());
            dual_faceverts.push_back(s2b.Center());
            } else {
            dual_faceverts.push_back(s2b.Center());
            dual_faceverts.push_back(s1.Center());
            }
            ++count;
        }
        m_dual_vertex_form_indices.push_back(count);
    }


    dual_faceindices.clear();
    dual_faceindices.resize(dual_faceverts.size());
    for(int i=0; i < dual_faceindices.size(); ++i) {
        dual_faceindices[i] = i;
    }




    std::cout << "Normalizing duals" << std::endl;

    mtao::normalizeToBBoxInPlace(dual_faceverts,m_bbox);
    mtao::normalizeToBBoxInPlace(dual_edgeverts,m_bbox);
    mtao::normalizeToBBoxInPlace(dual_vertices,m_bbox);







    //package->vertex_normals = m_mesh->getNormals();
    //package->vertex_normals = m_mesh->getNormals(NormalTriangleMesh::Area_Normal);
    //package->vertex_normals = m_mesh->getNormals(NormalTriangleMesh::Angle_Normal);
    package->vertex_normals = m_mesh->getNormals(NormalTriangleMesh::Mean_Curvature_Normal);






















    std::cout << "Done loading mesh metadata" << std::endl;
    emit meshLoaded(package);
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

        std::cout << "Got a key even!" << std::endl;
    switch(event->key()){
    default:
        m_glwidget->keyPressEvent(event);
        QMainWindow::keyPressEvent(event);
        return;
    }
    event->accept();
}
