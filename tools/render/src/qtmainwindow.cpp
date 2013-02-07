#include "../include/qtmainwindow.h"
#include <QMenu>
#include <QMenuBar>
#include <QAction>
#include <QFileDialog>
#include <QDir>

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
            , m_glwidget,SLOT(receiveMesh(const MeshPackage &)));
    setCentralWidget(m_glwidget);
}

void MainWindow::openFile(const QString & filename) {
    MeshPackage package;
    emit meshLoaded(package);
}

void MainWindow::openFile() {
    QFileDialog::Options options(QFileDialog::HideNameFilterDetails);
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Choose file or directory"),
                QDir::homePath());

    if (filename.isNull()) {
        return;
    }
    openFile(filename);
}
