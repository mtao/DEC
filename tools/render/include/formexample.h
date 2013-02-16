#ifndef EXAMPLEWIDGET_H
#define EXAMPLEWIDGET_H
#include "../../../include/dec.hpp"
#include "mainwindow.h"
#include "packages.h"
#include <QWidget>
class ExampleWidget: public MainWindow
{
    Q_OBJECT
    public:
    ExampleWidget(QWidget * parent=0);
    public slots:
    virtual void openFile(const QString & filename = "");
};
#endif
