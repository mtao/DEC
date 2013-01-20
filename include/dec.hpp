#ifndef DEC_H
#define DEC_H
#include "simplicialComplex.hpp"

enum FormType {PRIMAL, DUAL};
template <typename SimplicialComplex,FormType Type, int N>
class Form: public SimplicialComplex::NumTraits::DynamicVector
{
public:
    typedef typename SimplicialComplex::NumTraits NumTraits;
    typedef typename NumTraits::Scalar Scalar;
    typedef typename NumTraits::DynamicVector Parent;
    Form(): Parent(0) {}
    Form(const SimplicialComplex & sc)
        : Parent(
              (Type == PRIMAL) ?  sc.template numSimplices<N>() : sc.template numSimplices <SimplicialComplex::Dim-N>()
                                  )
    {
        static_assert(N <= SimplicialComplex::Dim,"Form can't be of higher dim than top dim of simplicial complex");
        this->setConstant(Scalar(0));
    }
    void init(const SimplicialComplex & sc)
    {
        resize((Type == PRIMAL) ?  sc.template numSimplices<N>() :sc.template numSimplices<SimplicialComplex::Dim-N>()
                                   );
        setConstant(Scalar(0));

    }
};


template <typename SimplicialComplex>
class FormFactory{
protected:
    FormFactory(const SimplicialComplex & sc): m_sc(sc) {}

    template <FormType Type, int N>
    Form<SimplicialComplex,Type,N> genForm()
    {
        return Form<SimplicialComplex,Type,N>(m_sc);
    }

private:
    const SimplicialComplex & m_sc;
};

template <typename SC, int TmD>
class HiddenOperatorContainer: public HiddenOperatorContainer<SC,TmD-1>
{

    protected:
    const static int TopD = SC::Dim;
    const static int D = TopD-TmD;
    typedef typename SC::NumTraits NumTraits;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
    typedef typename NumTraits::DynamicVector Vector;
    typedef HiddenOperatorContainer<SC,TmD-1> Parent;
    HiddenOperatorContainer(const SC & sc)
        : Parent(sc)
        , m_d(sc.template b<D+1>().transpose())
        , m_hodge(sc.template numSimplices<D>())
    {
        for(auto&& s: sc.template constSimplices<D>())
        {
            m_hodge(s.Index()) = s.DualVolume() / s.Volume();
        }

    }

    protected://Data
    SparseMatrix m_d;
    Vector m_hodge;
};

template <typename SC>
class HiddenOperatorContainer<SC,0>{
    protected:
    HiddenOperatorContainer(const SC &) {}
};

template <typename SC>
class OperatorContainer: public HiddenOperatorContainer<SC,SC::Dim> {
    protected:
    OperatorContainer(const SC & sc): HiddenOperatorContainer<SC,SC::Dim>(sc) {}

};

template <typename SC>
class DEC: public FormFactory<SC>, public OperatorContainer<SC>
{
public:
    typedef SC SimplicialComplex;
    typedef typename SC::NumTraits NumTraits;
    typedef typename NumTraits::SparseMatrix SparseMatrix;
    static const unsigned int Dim = SC::Dim;
    typedef FormFactory<SC> FF;
    typedef OperatorContainer<SC> OC;
    DEC(const SimplicialComplex & sc)
        : FF(sc)
        , OC(sc)
        , m_sc(sc)
    {
    }


    template <int N>
    const SparseMatrix & d()const{return HiddenOperatorContainer<SC,Dim-N>::m_d;}
    template <int N>
    const SparseMatrix & h()const{return HiddenOperatorContainer<SC,Dim-N>::m_hodge;}

    template <FormType Type, int N>
    const Form<SC,Type,N+1> d(const Form<SC,Type,N> & f)const{return Form<SC,Type,N+1>(d<N>()*f);}
    template <int N>
    const Form<SC,DUAL,Dim-N> h(const Form<SC,PRIMAL,N> & f)const{return f.cwiseProduct(h<N>());}
    template <int N>
    const Form<SC,PRIMAL,N> h(const Form<SC,DUAL,Dim-N> & f)const{return f.cwiseQuotient(h<N>());}


private:
    const SimplicialComplex & m_sc;
};

#endif
