#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/qtglwidget.h"
#include "../../../include/dec.hpp"
#include "packages.h"

class MainWindow: public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0);

protected:
    void keyPressEvent(QKeyEvent *);

public slots:
    void openFile();
    virtual void openFile(const QString & filename);
protected:
    GLWidget * m_glwidget;
    std::unique_ptr<TriangleMeshf> m_mesh;
    std::unique_ptr<DEC<TriangleMeshf,true> > m_dec;
    std::vector<unsigned int> m_dual_vertex_form_indices;
    template <typename Form>
    auto makeFormPackage(const QString & name, const Form & form)
    ->
        decltype(mtao::makeFormPackage(name, form, m_dual_vertex_form_indices))
    {
        return mtao::makeFormPackage(name, form, m_dual_vertex_form_indices);
    }
signals:
    void meshLoaded(std::shared_ptr<const MeshPackage> package);
    //void meshLoaded(const MeshPackage2 & package);
    void formLoaded(const FormPackage & package);
    void dataLoaded();
};

#endif
