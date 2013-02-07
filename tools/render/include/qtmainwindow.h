#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/glwidget.h"

class MainWindow: public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);

public slots:
    void openFile();
    void openFile(const QString & filename);
private:
    GLWidget * m_glwidget;
signals:
    void meshLoaded(const MeshPackage & package);
    void dataLoaded();
};

#endif
