
#ifndef TPZANALYTICSOLUTION
#define TPZANALYTICSOLUTION

#include "pzvec.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzfunction.h"

#include <string>


// typedef void (ExactFunc)(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);

struct TPZAnalyticSolution
{
    
    class TForce : public TPZFunction<STATE>
    {
        TPZAnalyticSolution *fAnalytic;
        
    public:
        TForce(TPZAnalyticSolution *root) : fAnalytic(root)
        {
        }
        
        /** @brief Simpler version of Execute method which does not compute function derivatives */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f){
            fAnalytic->Force(x,f);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder(){
            return 5;
        }

    };
    
    class TExactState : public TPZFunction<STATE>
    {
        TPZAnalyticSolution *fAnalytic;
        
    public:
        TExactState(TPZAnalyticSolution *root) : fAnalytic(root)
        {
            
        }
        
        /**
         * @brief Performs function computation
         * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
         * @param f function values
         * @param df function derivatives
         */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df)
        {
            fAnalytic->Solution(x,f,df);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder()
        {
            return 5;
        }

    };
    
    TPZAutoPointer<TPZFunction<STATE> > ForcingFunction()
    {
        return new TForce(this);
    };
    
    TPZAutoPointer<TPZFunction<STATE> > Exact()
    {
        return new TExactState(this);
    }
    
    virtual ~TPZAnalyticSolution()
    {
        
    }
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) = 0;
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) = 0;

};

#ifdef _AUTODIFF

struct TElasticityExample1 : public TPZAnalyticSolution
{
     enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EStretchy, EShear,Etest1,Etest2 };
    
     static EDefState fProblemType;
    

    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        TPZManVector<REAL,3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    virtual void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    
    virtual ~TElasticityExample1()
    {
        
    }
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad);

    template<class TVar>
    static void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);

    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);

};


struct TLaplaceExample1 : public TPZAnalyticSolution
{
    virtual ~TLaplaceExample1()
    {
        
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    

    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);
    
    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);

    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

};

struct TLaplaceExampleSmooth : public TPZAnalyticSolution
{
    virtual ~TLaplaceExampleSmooth()
    {
        
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);
    
    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);
    
    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZVec<TVar> &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }
    
};

#endif

#endif
