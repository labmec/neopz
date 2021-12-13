
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
    
    TPZAnalyticSolution() : fSignConvention(1){

    }

    TPZAnalyticSolution(const TPZAnalyticSolution &cp);
    
    TPZAnalyticSolution &operator=(const TPZAnalyticSolution &copy);

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
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f) override {
            fAnalytic->Force(x,f);
            for(auto &it:f) it *= fAnalytic->fSignConvention;
        }
        /** @brief Simpler version of Execute method which does not compute function derivatives */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &nada) override {
            fAnalytic->Force(x,f);
            for(auto &it:f) it *= fAnalytic->fSignConvention;
        }
        
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const override {
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
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) override
        {
            fAnalytic->Sigma(x,df);
            TPZFNMatrix<9,STATE> dsol(df);
            fAnalytic->Solution(x, f, dsol);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const override
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
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE> &df) override
        {
            fAnalytic->Solution(x,f,df);
        }
        /**
         * @brief Performs function computation
         * @param x point coordinate which is suppose to be in real coordinate system but can be in master coordinate system in derived classes.
         * @param f function values
         * @param df function derivatives
         */
        virtual void Execute(const TPZVec<REAL> &x, TPZVec<STATE> &f) override
        {
            TPZFNMatrix<9,STATE> df(3,3);
            fAnalytic->Solution(x,f,df);
        }
        /** @brief Polynomial order of this function. */
        /** In case of non-polynomial function it can be a reasonable approximation order. */
        virtual int PolynomialOrder() const override
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
    };
    
    std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)> ExactSolution() const
    {
        return [this](const TPZVec<REAL> &loc, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
        {
            this->Solution(loc,result,deriv);
        };
    }
    
    std::function<void (const TPZVec<REAL> &loc, TPZVec<STATE> &result)> ForceFunc() const
    {
        return [this](const TPZVec<REAL> &loc, TPZVec<STATE> &result)
        {
            this->Force(loc, result);
            for(auto &it:result) it *= this->fSignConvention;
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

    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const override
    {
        TPZManVector<STATE,3> xstate(3),locforce(2);
        for(int i=0; i<3; i++) xstate[i]=x[i];
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    virtual void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const;
    
    
    virtual ~TElasticity2DAnalytic()
    {
        
    }
    TElasticity2DAnalytic() : TPZAnalyticSolution(), fProblemType(ENone), fPlaneStress(1)
    {

    }

    TElasticity2DAnalytic(const TElasticity2DAnalytic &cp);

    TElasticity2DAnalytic &operator=(const TElasticity2DAnalytic &copy);

    static TPZAutoPointer<TPZFunction<STATE> > ConstitutiveLawFunction()
    {
        TPZAutoPointer<TPZFunction<STATE> > result;
        TPZDummyFunction<STATE> *dummy = new TPZDummyFunction<STATE>(TElasticity2DAnalytic::ElasticDummy,4);
        //dummy->SetPolynomialOrder(4);
        result = TPZAutoPointer<TPZFunction<STATE> >(dummy);
        return result;
    }

    void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const override;
    
    void Sigma(const TPZVec<Fad<STATE> > &x, TPZFMatrix<Fad<STATE> > &sigma) const;
    
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
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const override;
    
};

struct TElasticity3DAnalytic : public TPZAnalyticSolution
{
    enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1,Etest2, ETestShearMoment, ESphere };
    
    EDefState fProblemType = ENone;
    
    REAL fE = 1.;
    
    REAL fPoisson = 0.3;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const override
    {
        TPZManVector<STATE,3> locforce(3);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
        force[2] = -locforce[2];
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const override;
    
    TElasticity3DAnalytic() : TPZAnalyticSolution(), fProblemType(ENone)
    {

    }
    
    virtual ~TElasticity3DAnalytic()
    {
        
    }

    TElasticity3DAnalytic (const TElasticity3DAnalytic  &cp);

    TElasticity3DAnalytic  &operator=(const TElasticity3DAnalytic  &copy);

    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const override;

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
    
    enum EExactSol {ENone, EConst, EX, ESinSin, ECosCos,  EArcTan, EArcTanSingular,ESinDist, E10SinSin,E2SinSin, ESinSinDirNonHom,ESinMark,ESteklovNonConst,EGalvisNonConst,EBoundaryLayer,EBubble, EBubble2D,ESinCosCircle, EHarmonic, EHarmonic2,
    ESquareRootUpper, ESquareRootLower, ESquareRoot, ELaplace2D};
    
    int fDimension = 2;
    
    EExactSol fExact = EArcTan;
    
    TPZManVector<REAL,3> fCenter;

    TPZFNMatrix<9,REAL> fTensorPerm;

    TPZFNMatrix<9,REAL> fInvPerm;

    int fmaxIter = 15;

    static double gC;

    TLaplaceExample1() : fCenter(3,0.), TPZAnalyticSolution()
    {
        fTensorPerm = fInvPerm = {{1,0,0},{0,1,0},{0,0,1}};
    }

    TLaplaceExample1(TPZFNMatrix<9,REAL> K, TPZFNMatrix<9,REAL> invK): fCenter(3,0.), TPZAnalyticSolution()
    {
        fTensorPerm = K;
        fInvPerm =invK;
    }

    void setPermeabilyTensor(TPZFNMatrix<9,REAL> K, TPZFNMatrix<9,REAL> invK);

    virtual ~TLaplaceExample1()
    {
        fExact = ENone;
        fDimension = -1;
    }

    TLaplaceExample1(const TLaplaceExample1 &cp);

    TLaplaceExample1 &operator=(const TLaplaceExample1 &copy);
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const override;

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
    void DivSigma(const TPZVec<REAL> &x, TVar &divsigma) const;
    
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const override
    {
        STATE locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const override
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
    
    void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const override;
    
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const;
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const;
    
    template<class TVar>
    void Permeability(const TPZVec<TVar> &x, TVar &Elast) const;
    
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const override;
    
    template<class TVar>
    void DivSigma(const TPZVec<REAL> &x, TVar &divsigma) const;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const override
    {
        REAL locforce = 0;
        force[0] = locforce;
    }

};

struct TStokesAnalytic : public TPZAnalyticSolution
{
    
    enum MProblemType {EStokes, ENavierStokes, EOseen, ENavierStokesCDG, EOseenCDG, EBrinkman};
    
    enum EExactSol {ENone, ECavity,  EKovasznay, EKovasznayCDG, ESinCos, ENoFlow, ESinCos3D, EPconst, EObstacles, EOneCurve , ESinCosBDS, ESinCosBDS3D, EGatica3D, ECouplingSD, ECouplingNSD, EVugs2D, EVugs3D, EInfiltrationNS};
    
    int fDimension = 2;
    
    MProblemType fProblemType = EStokes;

    EExactSol fExactSol = ESinCos;
    
    REAL fvisco = 1.; //Viscosity

    REAL multRa = 1.; //No-Flow

    REAL Pi = M_PI;
        
    REAL fcBrinkman = 1.;
        
    TPZManVector<REAL,3> fCenter;
    
    TStokesAnalytic() : fCenter(3,0.), TPZAnalyticSolution()
    {
        
    }
    
    virtual ~TStokesAnalytic()
    {
        
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol) const override;
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &flux) const;

    template<class TVar>
    void pressure(const TPZVec<TVar> &x, TVar &p) const;
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad) const;

    template<class TVar>
    void Duxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &Du) const;
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const override;
    
    template<class TVar>
    void SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const;
    
    template<class TVar>
    void DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) const override;

    
};

#endif

