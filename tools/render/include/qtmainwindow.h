#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/qtglwidget.h"
#include "../../../include/dec.hpp"
#include "packages.h"
#include "formbar.h"


class MainWindow: public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0,FormBar * bar = 0);

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
    Eigen::AlignedBox<TriangleMeshf::Scalar, TriangleMeshf::EmbeddedDim> m_bbox;
    template <typename Form>
    auto makeFormPackage(const QString & name, const Form & form)
    ->
        decltype(mtao::makeFormPackage(name, form, m_dual_vertex_form_indices))
    {
        return mtao::makeFormPackage(name, form, m_dual_vertex_form_indices);
    }
    FormBar * m_formbar = 0;
signals:
    void meshLoaded(std::shared_ptr<const MeshPackage> package);
    void particlesLoaded(std::shared_ptr<VertexBufferObject> vbo);
    //void meshLoaded(const MeshPackage2 & package);
    void formLoaded(const FormPackage & package);
    void dataLoaded();
};

#endif
