#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/qtglwidget.h"

class MainWindow: public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);

protected:
    void keyPressEvent(QKeyEvent *);

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
