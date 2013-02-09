#include "../include/qtmainwindow.h"
#include "../../../include/dec.hpp"
#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QFileDialog>
#include <QDir>
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
    setCentralWidget(m_glwidget);
}

void MainWindow::openFile(const QString & filename) {
    MeshPackage package;
    //TODO: Memleak for now!!
    std::unique_ptr<TriangleMeshf> mesh(readOBJtoSimplicialComplex<float>(filename.toStdString()));
    DEC<TriangleMeshf> dec(*mesh);
    typedef typename TriangleMeshf::Vector Vector;
    const static int vecsize = sizeof(Vector)/sizeof(float);
    std::vector<Vector> verts = mesh->vertices();//copy here so we can normalize coordinates
    auto&& packed_indices = mesh->simplicesToArray();


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
    auto&& packed_edge_indices = mesh->simplicesToArray<1>();
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
