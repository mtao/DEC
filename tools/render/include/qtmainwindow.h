#ifndef MAINWINDOW_H
#define MAINWINDOW_H
#include <QMainWindow>
#include "../include/qtglwidget.h"
#include "../../../include/dec.hpp"

namespace mtao{
    template <typename Form>
        constexpr RenderType formToRendertype(const Form &) {
            return RenderType((Form::Traits::NOut == 0)
                *RT_VERT+
            (Form::Traits::NOut == 1)
                * RT_EDGE+
            (Form::Traits::NOut == 2)
                * RT_FACE);


        }
    template <typename Form>
        FormPackage makeFormPackage(const QString & name, const Form & form)
        {
            auto&& tmparr = formToRenderable(form);
            return  {name,formToRendertype(form), std::make_shared<VertexBufferObject>(
                    tmparr.data(),
                    (Form::Traits::NOut+1)*tmparr.size(),
                    GL_STATIC_DRAW,1)};
        }
};
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
signals:
    void meshLoaded(const MeshPackage * package);
    void formLoaded(const FormPackage & package);
    void dataLoaded();
};

#endif
