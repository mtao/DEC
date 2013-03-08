#ifndef LLT_H
#define LLT_H
#include <Eigen/Dense>


template <typename Matrix>
struct DenseLLT
{
    DenseLLT(const Matrix & A)
    {
            L=A.template triangularView<Eigen::Lower>();
            int i,j,k;
            for(i=0; i<A.rows(); ++i)
            {
                for(j=0; j<i; ++j)
                {
                    L(i,j) -= L.row(i).head(j).dot(L.row(j).head(j));
                    L(i,j)/=L(j,j);
                }
                L(i,i) -= L.row(i).head(i).squaredNorm();
                L(i,i) = sqrt(L(i,i));
            }

    }
    template <typename Vector>
    void solve(const Vector & b, Vector & x)
    {
        x=L.template triangularView<Eigen::Lower>().transpose().solve(
        L.template triangularView<Eigen::Lower>().solve(b)
                    );
    }
    private:
    Matrix L;
};



/*
#include "../linear.h"
template <typename Scalar, typename MatrixTag, unsigned int DIM>
class LLTEngine
{
    public:
        typedef SolverEngineTraits<LLTEngine<Scalar,MatrixTag,DIM> > Traits;
        friend Traits;
        typedef typename Traits::Matrix Matrix;
        typedef typename Traits::Vector Vector;
        LLTEngine(): tolerance(Traits::tolerance) {}
        LLTEngine(Scalar s): tolerance(s) {}
        void compute(Matrix & A)
        {
            L=A.template triangularView<Eigen::Lower>();
            int i,j,k;
            for(i=0; i<A.rows(); ++i)
            {
                //L(i,j) = \sqrt{A(i,j) - \sum_{0<=k<j} L(i,k)L(j,k)}
                for(j=0; j<i; ++j)//j<i
                {
                    for(k=0; k<j; ++k)
                    {
                        L(i,j)-=L(i,k)*L(j,k);
                    }
                    L(i,j)/=L(j,j);
                }
                //L(i,i) = A(i,i) - \sum_{0<=k<j} L(i,k)L(j,k)
                for(k=0; k<i; ++k)
                {
                    L(i,i)-=L(i,k)*L(i,k);
                }
                L(i,i) = sqrt(L(i,i));
            }

        }
        void solve(const Matrix & A, const Vector & b, Vector & x)
        {
            x.setZero();
            int k=0;
            int i,j;
            Vector dx;
            do
            {
                dx=A.diagonal().asDiagonal().inverse()*(b-A*x);
                x+=dx;
                ++k;
            } while (tolerance<  dx.norm());
        }
    protected:
        Matrix L;
        Scalar tolerance;
};

template <typename Scalar, unsigned int DIM>
class LLTEngine<Scalar,matrix_tags::sparse_tag,DIM>
{
    public:
        typedef SolverEngineTraits<LLTEngine<Scalar,matrix_tags::sparse_tag,DIM> > Traits;
        friend Traits;
        typedef typename Traits::Matrix Matrix;
        typedef typename Traits::Vector Vector;
        LLTEngine(): tolerance(Traits::tolerance) {}
        LLTEngine(Scalar s): tolerance(s) {}
        void compute(Matrix & A)
        {
            L=A.template triangularView<Eigen::Lower>();
            int i,j,k;
            for(i=0; i<A.rows(); ++i)
            {
                //L(i,j) = \sqrt{A(i,j) - \sum_{0<=k<j} L(i,k)L(j,k)}
                for(j=0; j<i; ++j)//j<i
                {
                    for(k=0; k<j; ++k)
                    {
                        L(i,j)-=L(i,k)*L(j,k);
                    }
                    L(i,j)/=L(j,j);
                }
                //L(i,i) = A(i,i) - \sum_{0<=k<j} L(i,k)L(j,k)
                for(k=0; k<i; ++k)
                {
                    L(i,i)-=L(i,k)*L(i,k);
                }
                L(i,i) = sqrt(L(i,i));
            }

        }
        void solve(const Matrix & A, const Vector & b, Vector & x)
        {
            x.setZero();
            int k=0;
            int i,j;
            Vector dx;
            do
            {
                dx=A.diagonal().asDiagonal().inverse()*(b-A*x);
                x+=dx;
                ++k;
            } while (tolerance<  dx.norm());
        }
    protected:
        Matrix L;
        Scalar tolerance;
};
*/
#endif
