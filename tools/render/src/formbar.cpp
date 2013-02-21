#include "../include/formbar.h"
#include <QGroupBox>
#include <QTreeWidget>
#include <QTreeWidgetItem>
#include <QVBoxLayout>

FormBar::FormBar(QWidget * parent): QWidget(parent) {
    setLayout(new QVBoxLayout);
    QGroupBox * groupbox = new QGroupBox(tr("Form Selector"));
    layout()->addWidget(groupbox);

    groupbox->setLayout(new QVBoxLayout);
    m_treewidget = new QTreeWidget(this);
    groupbox->layout()->addWidget(m_treewidget);

    m_treewidget->setColumnCount(3);
    QStringList headers;
    headers << "Title" << "Dimension" << "Primality";
    m_treewidget->setHeaderLabels(headers);
    connect(m_treewidget, SIGNAL(itemSelectionChanged(void)),
            this, SLOT(itemSelectionChanged(void)));
    

}

void FormBar::itemSelectionChanged() {
    emit clearForms();
    for(auto&& item: m_treewidget->selectedItems()) {
        emit enableForm(item->text(0));
    }
}

void FormBar::receiveForm(const FormPackage & package) {

    if(m_treewidget->findItems(package.title, Qt::MatchExactly, 0).size() != 0) {
        return;
    }



    QTreeWidgetItem * item = new QTreeWidgetItem(m_treewidget);
    item->setText(0,package.title);
    switch(package.type & ~RT_DUAL) {
        case RT_VERT:
            item->setText(1,tr("0"));
            break;
        case RT_EDGE:
            item->setText(1,tr("1"));
            break;
        case RT_FACE:
            item->setText(1,tr("2"));
            break;

    }
    item->setText(2,(package.type & RT_DUAL)?tr("DUAL"):tr("PRIMAL"));
    m_treewidget->insertTopLevelItem(0,item);
}
