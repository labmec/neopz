
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
    
    /// integer to correct for the sign convention of the forcing term
    int fSignConvention = 1;
    
    class TForce : public TPZFunction<STATE>
    {
        const TPZAnalyticSolution *fAnalytic;
        
    public:
        TForce(const TPZAnalyticSolution *root) : fAnalytic(root)
        {
        }
        
        virtual ~TForce()
        {
            
        }
        /** @brief Simpler version of Execute method which does not compute function derivatives */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f){
            fAnalytic->Force(x,f);
            for(auto &it:f) it *= fAnalytic->fSignConvention;
        }
        /** @brief Simpler version of Execute method which does not compute function derivatives */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &nada){
            fAnalytic->Force(x,f);
            for(auto &it:f) it *= fAnalytic->fSignConvention;
        }
        
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const {
            return 5;
        }

    };
    
    class Tensor : public TPZFunction<STATE>
    {
        TPZAnalyticSolution *fAnalytic;
        
    public:
        Tensor(TPZAnalyticSolution *root) : fAnalytic(root)
        {
            
        }
        
        virtual ~Tensor()
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
            fAnalytic->Sigma(x,df);
            TPZFNMatrix<9,STATE> dsol(df);
            fAnalytic->Solution(x, f, dsol);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const
        {
            return 5;
        }

    };
    
    class TExactState : public TPZFunction<STATE>
    {
        const TPZAnalyticSolution *fAnalytic;
        
    public:
        TExactState(const TPZAnalyticSolution *root) : fAnalytic(root)
        {
            
        }
        virtual ~TExactState()
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
        /**
         * @brief Performs function computation
         * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
         * @param f function values
         * @param df function derivatives
         */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f)
        {
            TPZFNMatrix<9,STATE> df(3,3);
            fAnalytic->Solution(x,f,df);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const
        {
            return 5;
        }
        
    };
    
    TPZAutoPointer<TPZFunction<STATE> > ForcingFunction() const
    {
        return new TForce(this);
    };
    
    TPZAutoPointer<TPZFunction<STATE> > Exact() const
    {
        return new TExactState(this);
    }
    
    std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> ExactSolution()
    {
        return [this](const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
        {
            this->Solution(loc,result,deriv);
        };
    }
    
    TPZAutoPointer<TPZFunction<STATE> > TensorFunction()
    {
        return new Tensor(this);
    }
    
    virtual ~TPZAnalyticSolution()
    {
        
    }
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const = 0;
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const = 0;

    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const = 0;
};

#ifdef _AUTODIFF

struct TElasticity2DAnalytic : public TPZAnalyticSolution
{
     enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1, Etest2, EThiago, EPoly,
         ESquareRootUpper, ESquareRootLower, ESquareRoot
     };
    
     EDefState fProblemType = EDispx;
    
    /// fPlaneStress = 1 -> plane stress
    /// fPlaneStress = 0 -> plane strain
    int fPlaneStress = 1;
    
    static REAL gE;
    
    static REAL gPoisson;
    
    static int gOscilatoryElasticity;

    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
    {
        TPZManVector<REAL,3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    virtual void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;
    
    
    virtual ~TElasticity2DAnalytic()
    {
        
    }
    
    static TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction()
    {
        TPZAutoPointer<TPZFunction<STATE> > result;
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(TElasticity2DAnalytic::ElasticDummy,4);
        //dummy->SetPolynomialOrder(4);
        result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
        return result;
    }

    void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const;
    
    void Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma) const;
    
    template<class TVar>
    void DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const;
    
    template<typename TVar1, typename TVar2>
    void uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &disp) const;
    
    template<typename TVar1, typename TVar2>
    void graduxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &grad) const {
        TPZManVector<Fad<TVar1>,3> xfad(x.size());
        for(int i=0; i<2; i++)
        {
            Fad<TVar1> temp = Fad<TVar1>(2,i,x[i]);
            xfad[i] = temp;

        }
        xfad[2] = x[2];
        TPZManVector<Fad<TVar2>,3> result(2);
        uxy(xfad,result);
        grad.Resize(2,2);
        for (int i=0; i<2; i++) {
            for (int j=0; j<2; j++)
            {
                grad(j,i) = result[i].d(j);
            }
        }
    }

    template<typename TVar>
    static void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);

    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;
    
};

struct TElasticity3DAnalytic : public TPZAnalyticSolution
{
    enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1,Etest2, ETestShearMoment, ESphere };
    
    EDefState fProblemType = ENone;
    
    REAL fE = 1.;
    
    REAL fPoisson = 0.3;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
    {
        TPZManVector<REAL,3> locforce(3);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
        force[2] = -locforce[2];
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;
    
    virtual ~TElasticity3DAnalytic()
    {
        
    }
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const;

    template<class TVar>
    void DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const;
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const;
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad) const;
    
    template<class TVar>
    void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu) const;
    
    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
};


struct TLaplaceExample1 : public TPZAnalyticSolution
{
    
    enum EExactSol {ENone, EConst, EX, ESinSin, ECosCos, EArcTan, EArcTanSingular,ESinDist, E10SinSin, ESinSinDirNonHom,ESinMark,ESteklovNonConst,EGalvisNonConst,EBoundaryLayer,EBubble};
    
    int fDimension = 2;
    
    EExactSol fExact = EArcTan;
    
    TPZManVector<REAL,3> fCenter;
    
    TLaplaceExample1() : fCenter(3,0.)
    {
        
    }
    
    virtual ~TLaplaceExample1()
    {
        
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;


    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const;
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const;
    
    template<class TVar>
    static void Permeability(const TPZVec<TVar> &x, TVar &Elast);

    static void PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    template<class TVar>
    void SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma) const;
    
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const
    {
        TPZManVector<STATE,3> xco(3);
        for (int i=0; i<3; i++) {
            xco[i] = x[i];
        }
        SigmaLoc<STATE>(xco,tensor);
    }

};

class TLaplaceExampleTimeDependent : public TPZAnalyticSolution
{
    
public:
    
    enum MProblemType {ENone, ELinear, ESin, ECos};
    
    MProblemType fProblemType = ESin;
    
    REAL fTime = 0.; // time
    
    REAL fDelt = 0.1; // timestep
    
    REAL fK = 1.; // permeability
    
    virtual ~TLaplaceExampleTimeDependent()
    {
        
    }
    
    void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;
    
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const;
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const;
    
    template<class TVar>
    void Permeability(const TPZVec<TVar> &x, TVar &Elast) const;
    
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const;
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma) const;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
    {
        REAL locforce = 0;
        force[0] = locforce;
    }

};

struct TStokes2DAnalytic : public TPZAnalyticSolution
{
    
    enum EExactSol {ENone, EStokes0, EStokesLimit1, EBrinkman1, EDarcyLimit1};
    
    int fDimension = 2;
    
    EExactSol fProblemType = ENone;
    
    REAL fvisco = 1.;
        
    TPZManVector<REAL,3> fCenter;
    
    TStokes2DAnalytic() : fCenter(3,0.)
    {
        
    }
    
    virtual ~TStokes2DAnalytic()
    {
        
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) const;
    
    template<typename TVar1, typename TVar2>
    void uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &flux) const;

    template<typename TVar1, typename TVar2>
    void pressure(const TPZVec<TVar1> &x, TVar2 &p) const;
    
    template<typename TVar1, typename TVar2>
    void graduxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &grad) const;

    template<typename TVar1, typename TVar2>
    void Duxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &Du) const;
    
    template<typename TVar1, typename TVar2>
    void SigmaLoc(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &sigma) const;
    
    template<typename TVar1, typename TVar2>
    void DivSigma(const TPZVec<TVar1> &x, TPZVec<TVar2> &divsigma) const;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const
    {
        TPZManVector<REAL,3> locforce(3);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
        force[2] = -locforce[2];
    }
    
    template<typename TVar1, typename TVar2>
    void Sigma(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const{
        TPZManVector<STATE,3> xco(3);
        for (int i=0; i<3; i++) {
            xco[i] = x[i];
        }
        SigmaLoc<STATE>(xco,tensor);
    }
    
    
    
};

#endif

#endif
