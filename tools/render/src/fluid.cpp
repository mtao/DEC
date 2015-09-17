#include "../../../include/io.hpp"
#include "../../../include/util.hpp"
#include "../../../include/render.hpp"
#include "../../../include/advection.hpp"
#include "../../../include/solvers/linear/conjugate_gradient/pcg.hpp"
#include "../../../tools/meshgeneration/sphere/include/sphere.hpp"
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
    void show();

protected:

    void initializeMesh();
    void keyPressEvent(QKeyEvent *);
private:
    typedef typename TriangleMeshf::Vector Vector;
    Vector v;
    void step(float dt=0.02);
    void pressure();
    decltype(m_dec->template genForm<FormType::Primal,1>()) velocity;
    decltype(m_dec->template genForm<FormType::Primal,2>()) m_pressure;
    decltype(m_dec->template genForm<FormType::Dual,0>()) m_pressure_dual;
    std::unique_ptr<std::vector<Particle<DECType> > > particles;
};

void FluidWidget::show() {
    MainWindow::show();
    emit loadingNewMesh();
    SphereMeshFactory<float> smf(4);
    m_mesh.reset(new MainWindow::MeshType(smf.faces(), smf.vertices()));
    initializeMesh();

}

FluidWidget::FluidWidget(QWidget * parent): MainWindow(parent){
    grabKeyboard();
}
void FluidWidget::initializeMesh() {
    MainWindow::initializeMesh();
    m_pressure = m_dec->template genForm<FormType::Primal,2>();
    velocity = m_dec->template genForm<FormType::Primal,1>();
    velocity.expr = decltype(velocity.expr)::Random(velocity.expr.rows());
    m_pressure.expr = decltype(m_pressure.expr)::Random(m_pressure.expr.rows());

    velocity.expr = decltype(velocity.expr)::Zero(velocity.expr.rows());
    velocity.expr(5) = 0.1;
    /*
    velocity.expr(50) = 1;
    velocity.expr(100) = 1;
    velocity.expr(150) = 1;
    */
    velocity.expr /= 1.0;
    //velocity.expr = decltype(velocity.expr)::Random(velocity.expr.rows());
    m_pressure.expr = decltype(m_pressure.expr)::Zero(m_pressure.expr.rows());
    emit formLoaded(makeFormPackage("p1", velocity));
    emit formLoaded(makeFormPackage("p2", m_pressure));



    particles.reset(new std::vector<Particle<DECType > >(/*m_mesh->numSimplices()*/1, Particle<DECType >(*m_dec, Vector::Random())));
    std::transform(m_mesh->simplices().begin(), m_mesh->simplices().begin()+1, particles->begin(), [&] (const decltype(m_mesh->simplices()[0]) & s)
            -> Particle<DECType > {
        return Particle<DECType >(*m_dec, s.Center(),&s);
    });

    std::vector<Vector> ps(particles->size());
    for(int i=0; i < ps.size(); ++i) {
        ps[i] = mtao::normalizeToBBox((*particles)[i].p(),m_bbox);
    }
    emit particlesLoaded(std::make_shared<VertexBufferObject>((void*)ps.data(),ps.size(),GL_STATIC_DRAW,3));




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
        step(0.02);
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
        velocitylines[2*i+1] = velocitylines[2*i]+vfield[i];
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

    Eigen::SparseMatrix<float> dhdh2 = m_dec->d(m_dec->h(m_dec->d(m_dec->template h<2>()))).expr;
    Eigen::SparseMatrix<float> hdh2 = m_dec->h(m_dec->d(m_dec->h<2>())).expr.eval();

    m_pressure= m_dec->d(velocity);

    Eigen::VectorXf rhs = m_pressure.expr;
    SparseCholeskyPreconditionedConjugateGradientSolve(dhdh2, rhs, m_pressure.expr);
    velocity.expr -= hdh2 * m_pressure.expr;

    /*

    Eigen::SparseMatrix<float> dhd2 =(m_dec->d(m_dec->h(m_dec->template d<0,FormType::Dual>()))).expr;

    m_pressure_dual = m_dec->h(m_dec->d(velocity));
    std::cout << dhd2.rows() << " " << dhd2.cols() << std::endl;

    Eigen::VectorXf rhs = m_pressure_dual.expr;
    SparseCholeskyPreconditionedConjugateGradientSolve(dhd2, rhs, m_pressure_dual.expr);

    m_pressure = m_dec->h(m_pressure_dual);
    //velocity.expr -= m_dec->h(m_dec->d(m_dec->h(m_pressure))).expr;
    velocity.expr -= (m_dec->d(m_pressure_dual)).expr;
    */

    emit formLoaded(makeFormPackage("rhs", m_pressure));



    //typename Eigen::ConjugateGradient<decltype(dhdh2), Eigen::Lower> chol;
    //typename Eigen::ConjugateGradient<decltype(dhdh2), Eigen::Lower, typename Eigen::SimplicialLDLT<decltype(dhdh2)> > chol;
    //chol.compute(dhdh2);
    //if(chol.info() != Eigen::Success) std::cout << "Failed at decomposition" << std::endl;
    //Eigen::VectorXf ret = chol.solve(m_pressure.expr);//solve poisson problem
    //if(chol.info() != Eigen::Success) std::cout << "Failed at solving" << std::endl;
    //std::cout << "CG Iterations: " << chol.iterations() << "/" << chol.maxIterations()<< " Error " << chol.error() << "/" << chol.tolerance()<< std::endl;

    //m_pressure.expr = ret;
    //velocity.expr -= dh2 * m_pressure.expr;
//    std::cout << "Error norm: " << (rhs - (dhdh2*m_pressure.expr)).norm() << "/(" << m_dec->d(velocity).expr.norm()<< ","<< m_pressure.expr.norm() << ")"<< std::endl;
    std::cout << "Final energy: " << velocity.expr.norm() << std::endl;
    std::cout << "Final Divergence: " << m_dec->d(velocity).expr.norm() << std::endl;

    //m_pressure.expr = m_dec->h<2,FormType::Dual>().expr * m_pressure.expr;
    m_pressure.expr *= 1000;
    emit formLoaded(makeFormPackage("m_pressure", m_pressure));


    emit formLoaded(makeFormPackage("projected velocity", velocity));


    std::vector<Vector> ps(particles->size());
    for(int i=0; i < ps.size(); ++i) {
        (*particles)[i].advectInPlace(velocity,0.02);
        ps[i] = mtao::normalizeToBBox((*particles)[i].p(),m_bbox);
    }
    emit particlesLoaded(std::make_shared<VertexBufferObject>((void*)ps.data(),ps.size(),GL_STATIC_DRAW,3));



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
