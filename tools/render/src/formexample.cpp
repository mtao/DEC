#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include "../../../include/render.hpp"
//#include "../../../include/advection.hpp"
#include "../include/qtmainwindow.h"
#include "../include/packages.h"
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <QKeyEvent>
class ExampleWidget: public MainWindow
{
public:
    ExampleWidget(QWidget * parent=0);
public:
    void openFile(const QString & filename);

protected:

    void keyPressEvent(QKeyEvent *);
private:
    typedef typename TriangleMeshf::Vector Vector;
    Vector v;
    void solveThings();
    decltype(m_dec->template genForm<PRIMAL_FORM,0>()) p0form;
    decltype(m_dec->template genForm<PRIMAL_FORM,1>()) p1form;
    decltype(m_dec->template genForm<PRIMAL_FORM,2>()) p2form;
    decltype(m_dec->template genForm<DUAL_FORM  ,0>()) d0form;
    decltype(m_dec->template genForm<DUAL_FORM  ,1>()) d1form;
    decltype(m_dec->template genForm<DUAL_FORM  ,2>()) d2form;
};
#include <iostream>
ExampleWidget::ExampleWidget(QWidget * parent): MainWindow(parent){
    grabKeyboard();
}
void ExampleWidget::openFile(const QString & filename) {
    MainWindow::openFile(filename);
    p2form = m_dec->template genForm<PRIMAL_FORM,2>();
    p1form = m_dec->template genForm<PRIMAL_FORM,1>();
    p0form = m_dec->template genForm<PRIMAL_FORM,0>();

    d2form = m_dec->template genForm<DUAL_FORM,2>();
    d1form = m_dec->template genForm<DUAL_FORM,1>();
    d0form = m_dec->template genForm<DUAL_FORM,0>();

    p0form.expr = decltype(p0form.expr)::Random(p0form.expr.rows());
    p1form.expr = decltype(p1form.expr)::Random(p1form.expr.rows());
    p2form.expr = decltype(p2form.expr)::Random(p2form.expr.rows());

    d0form.expr = decltype(d0form.expr)::Random(d0form.expr.rows());
    d1form.expr = decltype(d1form.expr)::Random(d1form.expr.rows());
    d2form.expr = decltype(d2form.expr)::Random(d2form.expr.rows());
    emit formLoaded(makeFormPackage("p0", p0form));
    emit formLoaded(makeFormPackage("p1", p1form));
    emit formLoaded(makeFormPackage("p2", p2form));
    emit formLoaded(makeFormPackage("d0", d0form));
    emit formLoaded(makeFormPackage("d1", d1form));
    emit formLoaded(makeFormPackage("d2", d2form));
    /*
    Particle<DEC<TriangleMeshf, true> > p(*m_dec,Vector(1,1,1));
    std::vector<Vector> ps;
    ps.resize(20);
    for(int i=1; i < ps.size(); ++i) {
        ps[i] = p.p();
    }
    //emit particlesLoaded(std::make_shared<VertexBufferObject>(&p.p()(0),1,GL_STATIC_DRAW,3));
    emit particlesLoaded(std::make_shared<VertexBufferObject>((void*)ps.data(),ps.size(),GL_STATIC_DRAW,3));
    */
}

void ExampleWidget::keyPressEvent(QKeyEvent * event) {

    switch(event->key()) {
    case Qt::Key_R:
        break;
    case Qt::Key_T:
        solveThings();
        break;
    }
}

#include <iostream>
void ExampleWidget::solveThings() {
    static std::default_random_engine generator;
    std::uniform_int_distribution<int> rand;
    p1form.expr = decltype(p2form.expr)::Zero(p1form.expr.rows());
    rand = std::uniform_int_distribution<int>(0,p1form.expr.rows()-1);
    p1form.expr(rand(generator)) = 1;
    emit formLoaded(makeFormPackage("initialVelocity", p1form));

    p2form = m_dec->d(p1form);
    emit formLoaded(makeFormPackage("rhs", p2form));

    Eigen::SparseMatrix<float> dhdh = m_dec->d(m_dec->h(m_dec->d(m_dec->template h<2>()))).expr.eval();

    typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::SimplicialLDLT<decltype(dhdh)> > chol;
    chol.setMaxIterations(dhdh.rows()/5);
    chol.setTolerance(0.001);
    //typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::CholmodBaseSupernodalLLT<decltype(dhdh)> > chol;
    chol.compute(dhdh);
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at dec->mposition" << std::endl;
    }
    Eigen::VectorXf ret = chol.solve(p2form.expr);//solve poisson problem
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at solving" << std::endl;
    }
    p2form.expr = ret / ret.lpNorm<Eigen::Infinity>();
    emit formLoaded(makeFormPackage("pressure", p2form));
    p1form.expr = m_dec->d(m_dec->h(p2form)).expr;
    p1form.expr /= p1form.expr.lpNorm<Eigen::Infinity>()*.2;
    emit formLoaded(makeFormPackage("projected velocity", p1form));

}


#ifdef FAKENESS

#define CHECK_INTERIORS

bool dumb = true;
void MainWindow::randomData() {
    if(!m_dec || !m_mesh) {return;}
    auto&& form2 = m_dec->template genForm<PRIMAL_FORM,2>();
    auto&& form1 = m_dec->template genForm<PRIMAL_FORM,1>();
    auto&& form0 = m_dec->template genForm<PRIMAL_FORM,0>();
#ifdef CHECK_INTERIORS
    auto&& interior2 = m_mesh->template interior<2>();
    for(int i=0; i < interior2.rows(); ++i) {

        form2.expr(i) = interior2.diagonal()(i);
        if(dumb)
            form2.expr(i) = 1-form2.expr(i);
    }
    auto&& interior1 = m_mesh->template interior<1>();
    for(int i=0; i < interior1.rows(); ++i) {
        form1.expr(i) = interior1.diagonal()(i);
        if(dumb)
            form1.expr(i) = 1-form1.expr(i);
    }
    auto&& interior0 = m_mesh->template interior<0>();
    for(int i=0; i < interior0.rows(); ++i) {
        form0.expr(i) = interior0.diagonal()(i);
        if(dumb)
            form0.expr(i) = 1-form0.expr(i);
    }
    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test0",form0));
    return;
#endif
    //form2.expr = decltype(form2.expr)::Random(form2.expr.rows());

    Eigen::SparseMatrix<float> dhdh = m_dec->h(m_dec->d(m_dec->h(m_dec->template d<1>()))).expr.eval();
    /*
    if(dumb) {
        form2 = m_2form;
        dumb = false;
    std::cout << dhdh << std::endl;
    std::cout << m_dec->template h<2>().expr << std::endl;
    std::cout << m_dec->template d<0,DUAL_FORM>().expr << std::endl;
    std::cout << m_dec->template h<1,DUAL_FORM>().expr << std::endl;
    std::cout << m_dec->template d<1,PRIMAL_FORM>().expr << std::endl;
    for(auto&& s: m_mesh->template simplices<1>()) {
        std::cout << s.Volume() << " " << s.DualVolume() << std::endl;
    }
    }
    */

    /*
    form2.expr = dhdh*form2.expr;
    std::cout << form2.expr.norm() << std::endl;
    form2.expr.normalize();

    */



    //typename Eigen::SimplicialLDLT<decltype(dhdh)> chol;
    typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::SimplicialLDLT<decltype(dhdh)> > chol;
    chol.setMaxIterations(dhdh.rows()/5);
    chol.setTolerance(0.1);
    //typename Eigen::ConjugateGradient<decltype(dhdh), Eigen::Lower, typename Eigen::CholmodBaseSupernodalLLT<decltype(dhdh)> > chol;
    chol.compute(dhdh);
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at dec->mposition" << std::endl;
    }
    Eigen::VectorXf ret = chol.solve(m_1form.expr);//solve poisson problem
    if(chol.info() != Eigen::Success)
    {
        std::cout << "Failed at solving" << std::endl;
    }
    m_1form.expr = ret;
    //std::cout << "Norm error: " << (dhdh* ret - m_1form.expr).norm() << std::endl;
    form1.expr = ret / ret.lpNorm<Eigen::Infinity>();
    /*
    static int WHICH = 0;
    form2.expr = m.col(WHICH);
    std::cout << form2.expr.transpose() << std::endl;
    WHICH = (WHICH+1)%100;
    */
    /*
    form2 = m_dec->d(m_dec->h(m_dec->d(m_dec->h(m_2form))));
    form2.expr = form2.expr / form2.expr.lpNorm<Eigen::Infinity>();
    */
    //form2 = m_dec->d(m_dec->h(form1));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));
    return;
    m_2form = form2;

    //form1 = m_dec->h(m_dec->d(m_dec->h(form2)));//apply codifferential operator
    //m_2form = form2;
    //form2 = m_2form;

    for(int i=0; i < form1.expr.rows(); ++i) {
        form1.expr(i) = m_mesh->template simplexByIndex<1>(i).DualVolume();
    }
    form1.expr = form1.expr / form1.expr.lpNorm<Eigen::Infinity>();

    m_glwidget->recieveForm(mtao::makeFormPackage("Test2",form2));
    m_glwidget->recieveForm(mtao::makeFormPackage("Test1",form1));

    for(int i=0; i < form0.expr.rows(); ++i)
    {
        form0.expr(i) = -1;
    }
    m_glwidget->recieveForm(mtao::makeFormPackage("Test0",form0));
}

void MainWindow::openFile() {
    QFileDialog::Options options(QFileDialog::HideNameFilterDetails);
    QString filename = QFileDialog::getOpenFileName(
                this,
                tr("Choose file or directory"),
                QDir::homePath(),
                tr("OBJ Files (*.obj)")
                );

    if (filename.isNull()) {
        return;
    }
    openFile(filename);
}
*/
#endif

#include <QApplication>
int main(int argc, char * argv[]) {
    QApplication a(argc,argv);
    ExampleWidget * mw = new ExampleWidget();
    QStringList args = a.arguments();
    if(args.size() > 1) {
        mw->openFile(args[1]);
    }
    mw->show();
    

    return a.exec();
}
