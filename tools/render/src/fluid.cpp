#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include "../../../include/render.hpp"
#include "../../../include/advection.hpp"
#include "../include/qtmainwindow.h"
#include "../include/packages.h"
#include <Eigen/SparseCholesky>
#include <Eigen/IterativeLinearSolvers>
#include <random>
#include <QKeyEvent>
#include <iostream>
class FluidWidget: public MainWindow
{
public:
    FluidWidget(QWidget * parent=0);
public:
    void openFile(const QString & filename);

protected:

    void keyPressEvent(QKeyEvent *);
private:
    typedef typename TriangleMeshf::Vector Vector;
    Vector v;
    void step(float dt=0.02);
    void pressure();
    decltype(m_dec->template genForm<PRIMAL_FORM,1>()) velocity;
    decltype(m_dec->template genForm<PRIMAL_FORM,2>()) m_pressure;
    std::unique_ptr<std::vector<Particle<DECType> > > particles;
};
FluidWidget::FluidWidget(QWidget * parent): MainWindow(parent){
    grabKeyboard();
}
void FluidWidget::openFile(const QString & filename) {
    MainWindow::openFile(filename);
    m_pressure = m_dec->template genForm<PRIMAL_FORM,2>();
    velocity = m_dec->template genForm<PRIMAL_FORM,1>();
    velocity.expr = decltype(velocity.expr)::Random(velocity.expr.rows());
    m_pressure.expr = decltype(m_pressure.expr)::Random(m_pressure.expr.rows());

    velocity.expr = decltype(velocity.expr)::Zero(velocity.expr.rows());
    velocity.expr(5) = 1;
    velocity.expr(50) = 1;
    velocity.expr(100) = 1;
    velocity.expr(150) = 1;
    //velocity.expr = decltype(velocity.expr)::Random(velocity.expr.rows());
    m_pressure.expr = decltype(m_pressure.expr)::Zero(m_pressure.expr.rows());
    emit formLoaded(makeFormPackage("p1", velocity));
    emit formLoaded(makeFormPackage("p2", m_pressure));



    /*
    particles.reset(new std::vector<Particle<DECType > >(m_mesh->numSimplices(), Particle<DECType >(*m_dec, Vector::Random())));
    std::transform(m_mesh->simplices().begin(), m_mesh->simplices().end(), particles->begin(), [&] (const decltype(m_mesh->simplices()[0]) & s)
            -> Particle<DECType > {
        return Particle<DECType >(*m_dec, s.Center(),&s);
    });
    std::vector<Vector> ps(particles->size());
    for(int i=0; i < ps.size(); ++i) {
        ps[i] = mtao::normalizeToBBox((*particles)[i].p(),m_bbox);
    }
    emit particlesLoaded(std::make_shared<VertexBufferObject>((void*)ps.data(),ps.size(),GL_STATIC_DRAW,3));
    */
    std::vector<Vector> vfield = m_dec->velocityField(velocity);
    m_glwidget->getVels(vfield);

    std::vector<Vector> velocitylines(2*vfield.size());
    for(int i=0; i < velocitylines.size()/2; ++i) {
        velocitylines[2*i] = mtao::normalizeToBBox(m_mesh->simplex(i).Center(),m_bbox);
        velocitylines[2*i+1] = velocitylines[2*i]+.01*vfield[i];
    }
    m_glwidget->getVels(velocitylines);
}

void FluidWidget::keyPressEvent(QKeyEvent * event) {

    switch(event->key()) {
    case Qt::Key_R:
        break;
    case Qt::Key_T:
        step(0.02);
        break;
    }
}

void FluidWidget::step(float dt) {
    emit formLoaded(makeFormPackage("initialVelocity", velocity));
    std::cout << "Initial energy: " << velocity.expr.norm() << std::endl;
    std::cout << "Initial Divergence: " << m_dec->d(velocity).expr.norm() << std::endl;
    pressure();
    //m_mesh->advect(velocity,.02);
    std::vector<Vector> vfield = m_dec->velocityField(velocity);
    m_glwidget->getVels(vfield);

    std::vector<Vector> velocitylines(2*vfield.size());
    for(int i=0; i < velocitylines.size()/2; ++i) {
        velocitylines[2*i] = mtao::normalizeToBBox(m_mesh->simplex(i).Center(),m_bbox);
        velocitylines[2*i+1] = velocitylines[2*i]+.01*vfield[i];
    }
    m_glwidget->getVels(velocitylines);
    /*

    std::vector<Vector> ps(particles->size());
    for(int i=0; i < ps.size(); ++i) {
        (*particles)[i].advectInPlace(velocity,0.02);
        ps[i] = mtao::normalizeToBBox((*particles)[i].p(),m_bbox);
    }
    emit particlesLoaded(std::make_shared<VertexBufferObject>((void*)ps.data(),ps.size(),GL_STATIC_DRAW,3));
    emit formLoaded(makeFormPackage("active", (*particles)[0].activeSimplex()));
    */
}
void FluidWidget::pressure() {

    emit formLoaded(makeFormPackage("rhs", m_pressure));
    Eigen::SparseMatrix<float> hdh2 = h(d(h<2>())).expr.eval();
    Eigen::SparseMatrix<float> hd2_ = m_dec->h(m_dec->h(m_dec->template d<2,DUAL_FORM>())).expr.eval();
    Eigen::SparseMatrix<float> d1 =m_dec->template d<1>().expr;
    Eigen::SparseMatrix<float> dhdh2 = d1*hdh2;
    m_pressure.expr = d1 * velocity.expr;

    //typename Eigen::ConjugateGradient<decltype(dhdh2), Eigen::Lower, typename Eigen::SimplicialLDLT<decltype(dhdh2)> > chol;
    typename Eigen::ConjugateGradient<decltype(dhdh2), Eigen::Lower> chol;
    chol.compute(dhdh2);
    if(chol.info() != Eigen::Success) std::cout << "Failed at decomposition" << std::endl;
    Eigen::VectorXf ret = chol.solve(m_pressure.expr);//solve poisson problem
    if(chol.info() != Eigen::Success) std::cout << "Failed at solving" << std::endl;
    std::cout << "CG Iterations: " << chol.iterations() << "/" << chol.maxIterations()<< " Error " << chol.error() << "/" << chol.tolerance()<< std::endl;

    m_pressure.expr = ret;
    velocity.expr -= hdh2 * m_pressure.expr;
    std::cout << "Error norm: " << (m_dec->d(velocity) - m_pressure).expr.norm() << "/(" << m_dec->d(velocity).expr.norm()<< ","<< m_pressure.expr.norm() << ")"<< std::endl;
    std::cout << "Final energy: " << velocity.expr.norm() << std::endl;
    std::cout << "Final Divergence: " << m_dec->d(velocity).expr.norm() << std::endl;


    m_pressure.expr = 1000*ret;
    emit formLoaded(makeFormPackage("m_pressure", m_pressure));


    emit formLoaded(makeFormPackage("projected velocity", velocity));





}


#include <QApplication>
int main(int argc, char * argv[]) {
    QApplication a(argc,argv);
    FluidWidget * mw = new FluidWidget();
    QStringList args = a.arguments();
    mw->show();

    if(args.size() > 1) {
        mw->openFile(args[1]);
    }

    return a.exec();
}
