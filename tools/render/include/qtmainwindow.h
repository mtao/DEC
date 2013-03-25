#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/qtglwidget.h"
#include "../../../include/dec.hpp"
#include "../../../include/trianglemesh.h"
#include "packages.h"
#include "formbar.h"


class MainWindow: public QMainWindow {
    Q_OBJECT
public:
    MainWindow(QWidget *parent = 0,FormBar * bar = 0);
    typedef NormalTriangleMesh MeshType;
    typedef DEC<MeshType, false > DECType;

protected:
    void keyPressEvent(QKeyEvent *);

public slots:
    void openFile();
    virtual void openFile(const QString & filename);
    virtual void initializeMesh();
protected:
    GLWidget * m_glwidget;
    std::unique_ptr<MeshType> m_mesh;
    std::unique_ptr<DECType> m_dec;
    std::vector<unsigned int> m_dual_vertex_form_indices;
    Eigen::AlignedBox<MeshType::Scalar, MeshType::EmbeddedDim> m_bbox;
    template <typename Form>
    auto makeFormPackage(const QString & name, const Form & form)
    ->
        decltype(mtao_internal::makeFormPackage(name, form, m_dual_vertex_form_indices))
    {
        return mtao_internal::makeFormPackage(name, form, m_dual_vertex_form_indices);
    }
    FormBar * m_formbar = 0;
signals:
    void meshLoaded(std::shared_ptr<const MeshPackage> package);
    void particlesLoaded(std::shared_ptr<VertexBufferObject> vbo);
    //void meshLoaded(const MeshPackage2 & package);
    void formLoaded(const FormPackage & package);
    void dataLoaded();
    void loadingNewMesh();
};

#endif
