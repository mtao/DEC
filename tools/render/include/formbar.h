#ifndef FORMBAR_H
#define FORMBAR_H

#include <QWidget>
#include "packages.h"
class QTreeWidget;

class FormBar: public QWidget{
    Q_OBJECT
    public:
        FormBar(QWidget * parent=0);
    protected:
        QTreeWidget * m_treewidget;
    public slots:
    void receiveForm(const FormPackage & package);
    void itemSelectionChanged();
signals:
    void enableForm(const QString & formname);
    void disableForm(const QString & name);
    void clearForms();

};


#endif
