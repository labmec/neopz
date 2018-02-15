
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
        
        virtual ~TForce()
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
        virtual int PolynomialOrder()
        {
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
    
    TPZAutoPointer<TPZFunction<STATE> > TensorFunction()
    {
        return new Tensor(this);
    }
    
    virtual ~TPZAnalyticSolution()
    {
        
    }
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force) = 0;
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) = 0;

    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) = 0;
};

#ifdef _AUTODIFF

struct TElasticity2DAnalytic : public TPZAnalyticSolution
{
     enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1, Etest2,
         ESquareRootUpper, ESquareRootLower, ESquareRoot
     };
    
     EDefState fProblemType = EDispx;
    
    /// fPlaneStress = 1 -> plane stress
    /// fPlaneStress = 0 -> plane strain
    int fPlaneStress = 1;
    
    REAL fE = 1.;
    
    REAL fPoisson = 0.3;

    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        TPZManVector<REAL,3> locforce(2);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
    }
    
    virtual void GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    
    virtual ~TElasticity2DAnalytic()
    {
        
    }
    
    void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma);
    
    void Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);
    
    template<typename TVar1, typename TVar2>
    void uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &disp) {

        if (fProblemType == Etest1) {
            disp[0] = TVar2(1. / 27.) * x[0] * x[0] * x[1] * x[1] * cos(TVar2(6. * M_PI) * x[0]) * sin(TVar2(7. * M_PI) * x[1]);
            disp[1] = TVar2(0.2) * exp(x[1]) * sin(TVar2(4. * M_PI) * x[0]);
        } else if (fProblemType == Etest2) {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] = ((TVar2(1) - x[0] * x[0])*(1. + x[1] * x[1] * x[1] * x[1]));
            disp[1] = ((TVar2(1) - x[1] * x[1])*(1. + x[0] * x[0] * x[0] * x[0]));
        }
        else if (fProblemType == ERot)//rotation
        {
            disp[0] = (TVar2) - x[1];
            disp[1] = (TVar2) x[0];
        }
        else if (fProblemType == EShear)//pure shear
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar2) x[1];
            disp[1] += (TVar2) 0.;
        } else if (fProblemType == EStretchx)//strech x
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar2) x[0];
            disp[1] += (TVar2) 0.;
        } else if (fProblemType == EUniAxialx) {
            if (fPlaneStress == 0) {
                disp[0] = x[0]*(1. - fPoisson * fPoisson) / fE;
                disp[1] = -x[1]*(1. + fPoisson) * fPoisson / fE;
            } else {
                disp[0] = x[0] / fE;
                disp[1] = -x[1] * fPoisson / fE;
            }
        } else if (fProblemType == EStretchy)//strech y
        {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar2) 0.;
            disp[1] += (TVar2) x[1];
        } else if (fProblemType == EDispx) {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += 1.;
            disp[0] += 0.;
        } else if (fProblemType == EDispy) {
            disp[0] = x[0]*0.;
            disp[1] = x[0]*0.;
            disp[0] += (TVar2) 0.;
            disp[0] += (TVar2) 1.;
        } else if (fProblemType == EBend) {
            TVar2 poiss = fPoisson;
            TVar2 elast = fE;
            if (fPlaneStress == 0) {
                poiss = poiss / (1. - poiss);
                elast /= (1 - fPoisson * fPoisson);
            }
            disp[0] = 5. * x[0] * x[1] / elast;
            disp[1] = (-poiss * 5. * x[1] * x[1] / 2. - 5. * x[0] * x[0] / 2.) / elast;
        } else if (fProblemType == ELoadedBeam) {
            TVar2 Est, nust, G;
            REAL MI = 5, h = 1.;
            G = fE / (2. * (1. + fPoisson));
            if (fPlaneStress == 0) {
                //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
                //            nust = fPoisson/(1.+fPoisson);
                Est = fE / ((1. - fPoisson * fPoisson));
                nust = fPoisson / (1 - fPoisson);
            } else {
                Est = fE;
                nust = fPoisson;
            }
            disp[0] = MI * h * h * x[1] / (2. * G) + MI * x[0] * x[0] * x[1] / (2. * Est) - MI * x[1] * x[1] * x[1] / (6. * G) + MI * nust * x[1] * x[1] * x[1] / (6. * Est);
            disp[1] = -MI * x[0] * x[0] * x[0] / (6. * Est) - MI * nust * x[0] * x[1] * x[1] / (2. * Est);
        } else if (fProblemType == ESquareRoot) {
#ifdef STATE_COMPLEX
            DebugStop();
#else
            TVar2 Est, nust, G, kappa;
            TVar2 theta = atan2(x[1], x[0]);
            TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
            G = fE / (2. * (1. + fPoisson));
            if (fPlaneStress == 0) {
                //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
                //            nust = fPoisson/(1.+fPoisson);
                Est = fE / ((1. - fPoisson * fPoisson));
                nust = fPoisson / (1 - fPoisson);
                kappa = 3. - 4. * fPoisson;
            } else {
                Est = fE;
                nust = fPoisson;
                kappa = (3. - fPoisson) / (1 + fPoisson);
            }
            TVar2 costh = cos(theta / 2.);
            TVar2 sinth = sin(theta / 2.);
            disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
            disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
            //        std::cout << "SQ x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
        } else if (fProblemType == ESquareRootLower) {
#ifdef STATE_COMPLEX
            DebugStop();
#else
            TVar2 Est, nust, G, kappa;
            TVar2 theta = atan2(x[1], x[0]);
            if (shapeFAD::val(theta) > 0.) {
                theta -= (2. * M_PI);
            }
            TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
            G = fE / (2. * (1. + fPoisson));
            if (fPlaneStress == 0) {
                //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
                //            nust = fPoisson/(1.+fPoisson);
                Est = fE / ((1. - fPoisson * fPoisson));
                nust = fPoisson / (1 - fPoisson);
                kappa = 3. - 4. * fPoisson;
            } else {
                Est = fE;
                nust = fPoisson;
                kappa = (3. - fPoisson) / (1 + fPoisson);
            }
            TVar2 costh = cos(theta / 2.);
            TVar2 sinth = sin(theta / 2.);
            disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
            disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
            //        std::cout << "SQL x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
        } else if (fProblemType == ESquareRootUpper) {
#ifdef STATE_COMPLEX
            DebugStop();
#else
            TVar2 Est, nust, G, kappa;
            TVar2 theta = atan2(x[1], x[0]);
            if (shapeFAD::val(theta) < 0.) {
                theta += (2. * M_PI);
            }
            TVar2 r = sqrt(x[0] * x[0] + x[1] * x[1]);
            G = fE / (2. * (1. + fPoisson));
            if (fPlaneStress == 0) {
                //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
                //            nust = fPoisson/(1.+fPoisson);
                Est = fE / ((1. - fPoisson * fPoisson));
                nust = fPoisson / (1 - fPoisson);
                kappa = 3. - 4. * fPoisson;
            } else {
                Est = fE;
                nust = fPoisson;
                kappa = (3. - fPoisson) / (1 + fPoisson);
            }
            TVar2 costh = cos(theta / 2.);
            TVar2 sinth = sin(theta / 2.);
            disp[0] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * costh * (kappa - 1. + 2. * sinth * sinth);
            disp[1] = 1 / (2. * G) * sqrt(r / (2. * M_PI)) * sinth * (kappa + 1. - 2. * costh * costh);
            //        std::cout << "SQU x " << x << " theta " << theta << " disp " << disp << std::endl;
#endif
        } else {
            DebugStop();
        }
    }
    
    template<typename TVar1, typename TVar2>
    void graduxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &grad){
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
    void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);

    static void ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv);
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
};

struct TElasticity3DAnalytic : public TPZAnalyticSolution
{
    enum EDefState  {ENone, EDispx, EDispy, ERot, EStretchx, EUniAxialx, EStretchy, EShear, EBend, ELoadedBeam, Etest1,Etest2, ETestShearMoment, ESphere };
    
    EDefState fProblemType = ENone;
    
    REAL fE = 1.;
    
    REAL fPoisson = 0.3;
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        TPZManVector<REAL,3> locforce(3);
        DivSigma(x, locforce);
        force[0] = -locforce[0];
        force[1] = -locforce[1];
        force[2] = -locforce[2];
    }
    
    virtual void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    virtual ~TElasticity3DAnalytic()
    {
        
    }
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor);

    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma);
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad);
    
    template<class TVar>
    void Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu);
    
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
    void SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce;
        DivSigma(x, locforce);
        force[0] = locforce;
    }

    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor)
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
    
    void Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu);
    
    
    template<class TVar>
    void uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp);
    
    template<class TVar>
    void graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad);
    
    template<class TVar>
    void Permeability(const TPZVec<TVar> &x, TVar &Elast);
    
    
    template<class TVar>
    void Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma);
    
    virtual void Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma);
    
    template<class TVar>
    void DivSigma(const TPZVec<TVar> &x, TVar &divsigma);
    
    virtual void Force(const TPZVec<REAL> &x, TPZVec<STATE> &force)
    {
        REAL locforce = 0;
        force[0] = locforce;
    }

};

#endif

#endif
