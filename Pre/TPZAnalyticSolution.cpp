
#include "TPZAnalyticSolution.h"

#include "pzgmesh.h"
#include "pzgengrid.h"
#include "pzgeoel.h"
#include "TPZRefPatternTools.h"
#include "pzcheckgeom.h"
#include "TPZVTKGeoMesh.h"

#include "pzanalysis.h"
#include "pzstepsolver.h"

#include "TPZMaterial.h"
#include "pzbuildmultiphysicsmesh.h"

#include "TPZSSpStructMatrix.h"

#ifndef USING_MKL
#include "pzskylstrmatrix.h"
#endif


#include "pzlog.h"

#ifdef _AUTODIFF
#include "fadType.h"


static Fad<REAL> atan2(Fad<REAL> y, Fad<REAL> x)
{
    int sz = x.size();
    Fad<REAL> result(sz,atan2(y.val(),x.val()));
    REAL r = x.val()*x.val()+y.val()*y.val();
    for (int i=0; i<sz; i++) {
        result.fastAccessDx(i) = (-y.val()*x.fastAccessDx(i) + x.val()*y.fastAccessDx(i))/r;
    }
    return result;
}

static FADFADREAL FADatan2(FADFADREAL y, FADFADREAL x)
{
    int sz = x.size();
    FADFADREAL result(sz,atan2(y.val(),x.val()));
    Fad<REAL> r = x.val()*x.val()+y.val()*y.val();
    for (int i=0; i<sz; i++) {
        result.fastAccessDx(i) = (-y.val()*x.fastAccessDx(i) + x.val()*y.fastAccessDx(i))/r;
    }
    return result;
}

static FADFADREAL FADsin(FADFADREAL x)
{
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    int sz = x.size();
    FADFADREAL sina(sz,sinaval);
    for (int i=0; i<sz; i++) {
        sina.fastAccessDx(i) = cosaval*x.dx(i);
    }
    return sina;
}

static FADFADREAL FADcos(FADFADREAL x)
{
    FADREAL_ sinaval = sin(x.val());
    FADREAL_ cosaval = cos(x.val());
    int sz = x.size();
    FADFADREAL cosa(sz,cosaval);
    for (int i=0; i<sz; i++) {
        cosa.fastAccessDx(i) = -sinaval*x.dx(i);
    }
    return cosa;
}

static FADFADREAL FADexp(FADFADREAL x)
{
    FADREAL_ expaval = exp(x.val());
    int sz = x.size();
    FADFADREAL expa(sz,expaval);
    for (int i=0; i<sz; i++) {
        expa.fastAccessDx(i) = expaval*x.dx(i);
    }
    return expa;
}

static FADFADREAL FADsqrt(FADFADREAL x)
{
    FADREAL_ fadres = sqrt(x.val());
    int sz = x.size();
    FADFADREAL resa(sz,fadres);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = REAL(0.5)/fadres*x.dx(i);
    }
    return resa;
}

static FADFADREAL FADatan(FADFADREAL x)
{
    FADREAL_ fadres = atan(x.val());
    int sz = x.size();
    FADFADREAL resa(sz,fadres);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = 1./(1+x.val()*x.val())*x.dx(i);
    }
    return resa;
}

static FADFADREAL FADcosh(FADFADREAL x)
{
    FADREAL_ coshval = cosh(x.val());
    FADREAL_ sinhval = sinh(x.val());
    int sz = x.size();
    FADFADREAL resa(sz,coshval);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = sinhval*x.dx(i);
    }
    return resa;
}

static FADFADREAL FADsinh(FADFADREAL x)
{
    FADREAL_ coshval = cosh(x.val());
    FADREAL_ sinhval = sinh(x.val());
    int sz = x.size();
    FADFADREAL resa(sz,sinhval);
    for (int i=0; i<sz; i++) {
        resa.fastAccessDx(i) = coshval*x.dx(i);
    }
    return resa;
}

static const REAL FI = 1.;
static const REAL a = 0.5;
static const REAL b = 0.5;

REAL TElasticity2DAnalytic::gE = 1.;

REAL TElasticity2DAnalytic::gPoisson = 0.3;

int TElasticity2DAnalytic::gOscilatoryElasticity = 0;

template<>
void TElasticity2DAnalytic::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) const
{
    typedef FADFADREAL TVar;
    if(fProblemType == Etest1)
    {   
        FADFADREAL tmp = (FADFADREAL)(1./27.)*x[0]*x[0]*x[1]*x[1];
        disp[0] = tmp*FADcos((FADFADREAL)(6.*M_PI)*x[0])*FADsin((FADFADREAL)(7.*M_PI)*x[1]);
        disp[1] = (FADFADREAL)(0.2)*FADexp(x[1])*FADsin((FADFADREAL)(4.*M_PI)*x[0]);
    
    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1-x[0]*x[0])*(1+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1-x[1]*x[1])*(1+x[0]*x[0]*x[0]*x[0]));
    }
      
    else if(fProblemType ==ERot)//rotation
    {      
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] =(FADFADREAL)-x[1];
        disp[1] =(FADFADREAL) x[0];
      
    }
    
    else if(fProblemType == EShear)//pure shear
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) x[1];
        disp[1] += (FADFADREAL) 0. ;
    }
    else if(fProblemType == EStretchx)//strech x
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) x[0];
        disp[1] += (FADFADREAL) 0.;
    }
    else if(fProblemType == EUniAxialx)
    {
        if (fPlaneStress == 0) {
            disp[0] = x[0]*(1.-gPoisson*gPoisson)/gE;
            disp[1] = -x[1]*(1.+gPoisson)*gPoisson/gE;
        }
        else
        {
            disp[0] = x[0]/gE;
            disp[1] = -x[1]*gPoisson/gE;
        }
    }
    else if(fProblemType ==EStretchy)//strech y
    {    
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) 0.;
        disp[1] += (FADFADREAL) x[1];    
    }
    else if(fProblemType==EDispx)
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[0] +=   0.;
    }
    else if(fProblemType==EDispy)
    {     
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) 0.;
        disp[0] += (FADFADREAL) 1.;
    }
    else if (fProblemType == EThiago){
        disp[0] = FADcos(M_PI * x[0]) * FADsin(2 * M_PI * x[1]);
        disp[1] = FADcos(M_PI * x[1]) * FADsin(M_PI * x[0]);
    } else if(fProblemType == EPoly){
        disp[0] = 0. * x[0]; //*x[0]*x[1];
        disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
    } else if(fProblemType==EBend)
    {
        typedef TVar FADFADREAL;
        TVar poiss = gPoisson;
        TVar elast = gE;
        if(fPlaneStress == 0)
        {
            poiss = poiss/(1.-poiss);
            elast /= (1-gPoisson*gPoisson);
        }
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est,nust,G;
        REAL MI = 5, h = 1.;
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1-gPoisson);
//            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
//            nust = gPoisson/(1.+gPoisson);
        }
        else
        {
            Est = gE;
            nust = gPoisson;
        }
        disp[0] = MI*h*h*x[1]/(2.*G)+ MI*x[0]*x[0]*x[1]/(2.*Est)-MI *x[1]*x[1]*x[1]/(6.*G)+MI*nust*x[1]*x[1]*x[1]/(6.*Est);
        disp[1] = -MI*x[0]*x[0]*x[0]/(6*Est)-MI*nust*x[0]*x[1]*x[1]/(2.*Est);
    }
    else if(fProblemType == ESquareRoot)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    else if(fProblemType == ESquareRootLower)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        if (shapeFAD::val(theta) > 0.) {
            theta -= (2.*M_PI);
        }
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    else if(fProblemType == ESquareRootUpper)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        if (shapeFAD::val(theta) < 0.) {
            theta += (2.*M_PI);
        }
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = gE/(2.*(1.+gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*fE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE/((1.-gPoisson*gPoisson));
            nust = gPoisson/(1-gPoisson);
            kappa = 3.-4.*gPoisson;
        }
        else
        {
            Est = gE;
            nust = gPoisson;
            kappa = (3.-gPoisson)/(1+gPoisson);
        }
        TVar costh = FADcos(theta/2.);
        TVar sinth = FADsin(theta/2.);
        disp[0] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*costh*(kappa-1.+2.*sinth*sinth);
        disp[1] = 1/(2.*G)*FADsqrt(r/(2.*M_PI))*sinth*(kappa+1.-2.*costh*costh);
    }
    else
    {
        DebugStop();
    }
}

template<typename TVar1, typename TVar2>
void TElasticity2DAnalytic::uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &disp) const {
    
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
            disp[0] = x[0]*(1. - gPoisson * gPoisson) / gE;
            disp[1] = -x[1]*(1. + gPoisson) * gPoisson / gE;
        } else {
            disp[0] = x[0] / gE;
            disp[1] = -x[1] * gPoisson / gE;
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
    } else if (fProblemType == EThiago){
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        //disp[0] = ((1-x[0]*x[0])*(1+x[1]*x[1]*x[1]*x[1]));
        //disp[1] = ((1-x[1]*x[1])*(1+x[0]*x[0]*x[0]*x[0]));
        disp[0] = (TVar2) cos(M_PI * x[0])*(TVar2) sin(2 * M_PI * x[1]);
        disp[1] = (TVar2) cos(M_PI * x[1])*(TVar2) sin(M_PI * x[0]);
    } else if(fProblemType == EPoly){
        disp[0] = 0. * x[0]; //*x[0]*x[1];
        disp[1] = x[1] * x[0]; //(x[1]-1.)*(x[1]-x[0])*(x[0]+4.*x[1]);
    } else if (fProblemType == EBend) {
        TVar2 poiss = gPoisson;
        TVar2 elast = gE;
        if (fPlaneStress == 0) {
            poiss = poiss / (1. - poiss);
            elast /= (1 - gPoisson * gPoisson);
        }
        disp[0] = 5. * x[0] * x[1] / elast;
        disp[1] = (-poiss * 5. * x[1] * x[1] / 2. - 5. * x[0] * x[0] / 2.) / elast;
    } else if (fProblemType == ELoadedBeam) {
        TVar2 Est, nust, G;
        REAL MI = 5, h = 1.;
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
        } else {
            Est = gE;
            nust = gPoisson;
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
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
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
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
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
        G = gE / (2. * (1. + gPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*gPoisson)/((1+gPoisson)*(1.+gPoisson))*gE;
            //            nust = gPoisson/(1.+gPoisson);
            Est = gE / ((1. - gPoisson * gPoisson));
            nust = gPoisson / (1 - gPoisson);
            kappa = 3. - 4. * gPoisson;
        } else {
            Est = gE;
            nust = gPoisson;
            kappa = (3. - gPoisson) / (1 + gPoisson);
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

template
void TElasticity2DAnalytic::uxy<REAL,REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &disp) const;

#if (REAL != STATE)

template
void TElasticity2DAnalytic::uxy<REAL,STATE>(const TPZVec<REAL> &x, TPZVec<STATE> &disp) const;

template
void TElasticity2DAnalytic::uxy<STATE,STATE>(const TPZVec<STATE> &x, TPZVec<STATE> &disp) const;

#endif

template<class TVar>
void TElasticity2DAnalytic::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu)
{
    if(gOscilatoryElasticity != 0)
    {
        Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    }
    else
    {
        Elast = gE;
    }
    nu = TVar(gPoisson);
}

//template<>
//void TElasticity2DAnalytic::Elastic(const TPZVec<double> &x, double &Elast, double &nu)
//{
//  Elast = 1000.;
//    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
//    Elast = 1.;
//    nu = 0.;
//}

void TElasticity2DAnalytic::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    //TPZManVector<STATE> xstate(x.size());
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    
    STATE E = 1,nu = 0.3;
    Elastic(xstate,E,nu);
    result[0] = E;
    result[1] = nu;
}


template<>
void TElasticity2DAnalytic::graduxy(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &grad) const
{
    TPZManVector<Fad<Fad<REAL> >,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<Fad<REAL> > temp = Fad<Fad<REAL> >(3,Fad<REAL>(3,0.));
//      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        temp.val()= x[i];
        Fad<REAL> temp3(3,0.);
        for(int j=0; j<3; j++)
        {
            temp.fastAccessDx(j) = temp3;
        }
        Fad<REAL> temp2(3,1.);
        temp.fastAccessDx(i) = temp2;
//      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;    
//      xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<REAL> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(j,i) = result[i].d(j);
        }
    }
}

template
void TElasticity2DAnalytic::graduxy<REAL,REAL>(const TPZVec<REAL> &x, TPZFMatrix<REAL> &grad) const;

#if (REAL != STATE)

template
void TElasticity2DAnalytic::graduxy<REAL,STATE>(const TPZVec<REAL> &x, TPZFMatrix<STATE> &grad) const;

template
void TElasticity2DAnalytic::graduxy<STATE,STATE>(const TPZVec<STATE> &x, TPZFMatrix<STATE> &grad) const;

#endif

void TElasticity2DAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const {
    TPZManVector<STATE> xst(3);
    for(int i=0; i<3; i++) xst[i] = x[i];
    uxy(xst,u);
    graduxy(xst,gradu);
}

void TElasticity2DAnalytic::GradU(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,2);
    u[0] = result[0].val();
    u[1] = result[1].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

void TElasticity2DAnalytic::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
{
    TPZFNMatrix<4,STATE> grad;
    REAL E, nu;
    sigma.Resize(2,2);
    Elastic(x, E, nu);
    graduxy(x,grad);
    //uxy(x,u);
    if (fPlaneStress == 0)
    {
        STATE Fac = E/((STATE)1.+nu)/((STATE(1.)-STATE(2.)*nu));
        sigma(0,0) = Fac*((STATE(1.)-nu)*grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*((STATE(1.)-nu)*grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(STATE(2.)*(STATE(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
    else
    {
        STATE Fac = E/((STATE)1.-nu*nu);
        sigma(0,0) = Fac*(grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*(grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(STATE(2.)*(STATE(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
}

void TElasticity2DAnalytic::Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma) const {
    TPZFNMatrix<4,Fad<REAL> > grad;
    sigma.Resize(2,2);
    Fad<REAL>  E, nu;
    Elastic(x, E, nu);
    graduxy(x,grad);
    if (fPlaneStress == 0)
    {
        Fad<REAL>  Fac = E/(Fad<REAL>(1.)+nu)/((Fad<REAL>(1.)-Fad<REAL>(2.)*nu));
        sigma(0,0) = Fac*((Fad<REAL>(1.)-nu)*grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*((Fad<REAL>(1.)-nu)*grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(Fad<REAL>(2.)*(Fad<REAL>(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);
    }
    else
    {
        typedef Fad<REAL> TVar;
        TVar Fac = E/((TVar)1.-nu*nu);
        sigma(0,0) = Fac*(grad(0,0)+nu*grad(1,1));
        sigma(1,1) = Fac*(grad(1,1)+nu*grad(0,0));
        sigma(0,1) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,1)+grad(1,0));
        sigma(1,0) = sigma(0,1);

    }
}

template<class TVar>
//void TElasticity2DAnalytic::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma) const
void TElasticity2DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<4, Fad<TVar> > sigma(2,2);
    Sigma(xfad,sigma);
    divsigma[0] = sigma(0,0).dx(0)+sigma(0,1).dx(1);
    divsigma[1] = sigma(1,0).dx(0)+sigma(1,1).dx(1);
}



template
void TElasticity2DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<REAL> &divsigma) const;



template<class TVar>
void TElasticity3DAnalytic::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    
    if(fProblemType == Etest1)
    {
        disp[0] = TVar(1./27.)*x[0]*x[0]*x[1]*x[1]*cos(TVar(6.*M_PI)*x[0])*sin(TVar(7.*M_PI)*x[1]);
        disp[1] = TVar(0.2)*exp(x[1])*sin(TVar(4.*M_PI)*x[0]);
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1.-x[0]*x[0])*(1.+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1.-x[1]*x[1])*(1.+x[0]*x[0]*x[0]*x[0]));
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType ==ERot)//rotation
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = (TVar)-x[1];
        disp[1] = (TVar)x[0];
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType == EShear)//pure shear
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) x[1];
        disp[1] += (TVar) 0. ;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EStretchx)//strech x
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) x[0];
        disp[1] += (TVar) 0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EUniAxialx)//strech x
    {
        disp[0] += (TVar) x[0]/fE;
        disp[1] += (TVar) -x[1]*fPoisson/fE;
        disp[2] = -x[2]*fPoisson/fE;
    }
    else if(fProblemType ==EStretchy)//strech y
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) 0.;
        disp[1] += (TVar) x[1];
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispx)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[0] +=   0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispy)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (TVar) 0.;
        disp[0] += (TVar) 1.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EBend)
    {
        TVar poiss = fPoisson;
        TVar elast = fE;
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        disp[2] = x[0]*0.;
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est = fE,nust = fPoisson,G;
        REAL MI = 5, h = 1.;
        G = fE/(2.*(1.+fPoisson));
        // returning the solution of plane strain
        Est = fE/((1.-fPoisson*fPoisson));
        nust = fPoisson/(1-fPoisson);
        int offset = 2;
        int i0 = offset;
        int i1 = (offset+1)%3;
        int i2 = (offset+2)%3;
        
        disp[i0] = MI*h*h*x[i1]/(2.*G)+MI * x[i0]*x[i0]*x[i1]/(2.*Est)-MI *x[i1]*x[i1]*x[i1]/(6.*G)+MI*nust*x[i1]*x[i1]*x[i1]/(6.*Est);
        disp[i1] = -MI*x[i0]*x[i0]*x[i0]/(6.*Est)-MI*nust*x[i0]*x[i1]*x[i1]/(2.*Est);
        disp[i2] = x[i2]*0.;
    }
    else if(fProblemType == ETestShearMoment)
    {
        disp[0] = -((TVar) FI/fE*fPoisson) *x[0]*x[1]*x[2];
        disp[1] = ((TVar) FI/fE) *((TVar)(fPoisson/2.)*(x[0]*x[0]-x[1]*x[1])*x[2]-((TVar)(1./6.)*x[2]*x[2]*x[2]));
        TVar disp2A = FI/fE*(1./2.* x[1]*(fPoisson*x[0]*x[0]+x[2]*x[2])+1./6.*fPoisson*x[1]*x[1]*x[1]+(1.+fPoisson)*(b*b*x[1]-1./3.*x[1]*x[1]*x[1])-1./3.*a*a*fPoisson*x[1]);
        TVar series = (TVar) 0.;
        TVar minusone = (TVar) -1;
        for (int i=1; i<5; i++) {
            series += (minusone/TVar(i*i*i)*cos(i*M_PI*x[0]/a)*sinh(i*M_PI*x[1]/a)/cosh(i*M_PI*b/a));
            minusone *= (TVar) -1.;
        }
        series *= (-4.*a*a*a*fPoisson/(M_PI*M_PI*M_PI));
        disp[2] = disp2A+series;
    }
    else if(fProblemType == ESphere)
    {
        TPZManVector<TVar,3> xc(x);
        TVar radius2 = xc[0]*xc[0]+xc[1]*xc[1]+xc[2]*xc[2];
        TVar radius = sqrt(radius2);
        TVar acube = 10.*10.*10.;
        TVar bcube = 50.*50.*50.;
        TVar pi = 6.;
        TVar ur = pi*acube/(2.*fE*(bcube-acube)*radius2)*(2.*(1.-2.*fPoisson)*radius2*radius+(1.+fPoisson)*bcube);
        TVar cosphi = xc[2]/radius;
        TVar xyradius = sqrt(xc[0]*xc[0]+xc[1]*xc[1]);
        TVar sinphi = xyradius/radius;
        TVar costheta = xc[0]/xyradius;
        TVar sintheta = xc[1]/xyradius;
        disp[2] = ur*cosphi;
        disp[0] = ur*sinphi*costheta;
        disp[1] = ur*sinphi*sintheta;
    }
    else{
        DebugStop();
    }
}


template<>
void TElasticity3DAnalytic::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) const
{
    typedef FADFADREAL TVar;
    if(fProblemType == Etest1)
    {
        FADFADREAL tmp = (FADFADREAL)(1./27.)*x[0]*x[0]*x[1]*x[1];
        disp[0] = tmp*FADcos((FADFADREAL)(6.*M_PI)*x[0])*FADsin((FADFADREAL)(7.*M_PI)*x[1]);
        disp[1] = (FADFADREAL)(0.2)*FADexp(x[1])*FADsin((FADFADREAL)(4.*M_PI)*x[0]);
        disp[2] = x[0]*TVar(0.);

    }
    else if(fProblemType == Etest2)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] = ((1-x[0]*x[0])*(1+x[1]*x[1]*x[1]*x[1]));
        disp[1] = ((1-x[1]*x[1])*(1+x[0]*x[0]*x[0]*x[0]));
        disp[2] = x[0]*TVar(0.);
    }
    
    else if(fProblemType ==ERot)//rotation
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] =(FADFADREAL)-x[1];
        disp[1] =(FADFADREAL) x[0];
        disp[2] = x[0]*TVar(0.);

    }
    
    else if(fProblemType == EShear)//pure shear
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) x[1];
        disp[1] += (FADFADREAL) 0. ;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EStretchx)//strech x
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) x[0];
        disp[1] += (FADFADREAL) 0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType == EUniAxialx)//strech x
    {
        disp[0] += (TVar) x[0]/fE;
        disp[1] += (TVar) -x[1]*fPoisson/fE;
        disp[2] = -x[2]*fPoisson/fE;
    }
    else if(fProblemType ==EStretchy)//strech y
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) 0.;
        disp[1] += (FADFADREAL) x[1];
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispx)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] +=   1.;
        disp[0] +=   0.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EDispy)
    {
        disp[0] = x[0]*0.;
        disp[1] = x[0]*0.;
        disp[0] += (FADFADREAL) 0.;
        disp[0] += (FADFADREAL) 1.;
        disp[2] = x[0]*TVar(0.);
    }
    else if(fProblemType==EBend)
    {
        TVar poiss = fPoisson;
        TVar elast = fE;
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        disp[2] = x[0]*0.;
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est = fE,nust = fPoisson,G;
        REAL MI = 5, h = 1.;
        G = fE/(2.*(1.+fPoisson));
        // returning the solution of plane strain
        Est = fE/((1.-fPoisson*fPoisson));
        nust = fPoisson/(1-fPoisson);
        int offset = 2;
        int i0 = offset;
        int i1 = (offset+1)%3;
        int i2 = (offset+2)%3;

        disp[i0] = MI*h*h*x[i1]/(2.*G)+MI * x[i0]*x[i0]*x[i1]/(2.*Est)-MI *x[i1]*x[i1]*x[i1]/(6.*G)+MI*nust*x[i1]*x[i1]*x[i1]/(6.*Est);
        disp[i1] = -MI*x[i0]*x[i0]*x[i0]/(6*Est)-MI*nust*x[i0]*x[i1]*x[i1]/(2.*Est);
        disp[i2] = x[i2]*0.;
    }
    else if(fProblemType == ETestShearMoment)
    {
        disp[0] = -((TVar) FI/fE*fPoisson) *x[0]*x[1]*x[2];
        disp[1] = ((TVar) FI/fE) *((TVar)(fPoisson/2.)*(x[0]*x[0]-x[1]*x[1])*x[2]-((TVar)(1./6.)*x[2]*x[2]*x[2]));
        TVar disp2A = FI/fE*(1./2.* x[1]*(fPoisson*x[0]*x[0]+x[2]*x[2])+1./6.*fPoisson*x[1]*x[1]*x[1]+(1.+fPoisson)*(b*b*x[1]-1./3.*x[1]*x[1]*x[1])-1./3.*a*a*fPoisson*x[1]);
        TVar series = (TVar) 0.;
        TVar minusone = (TVar) -1;
        for (int i=1; i<5; i++) {
            series += (minusone/(i*i*i)*FADcos(i*M_PI*x[0]/a)*FADsinh(i*M_PI*x[1]/a)/FADcosh(i*M_PI*b/a));
            minusone *= (TVar) -1.;
        }
        series *= (-4.*a*a*a*fPoisson/(M_PI*M_PI*M_PI));
        disp[2] = disp2A+series;
    }
    else if(fProblemType == ESphere)
    {
        TPZManVector<TVar,3> xc(x);
        TVar radius2 = xc[0]*xc[0]+xc[1]*xc[1]+xc[2]*xc[2];
        TVar radius = FADsqrt(radius2);
        TVar acube = 10.*10.*10.;
        TVar bcube = 50.*50.*50.;
        TVar pi = 6.;
        TVar ur = pi*acube/(2.*fE*(bcube-acube)*radius2)*(2.*(1.-2.*fPoisson)*radius2*radius+(1.+fPoisson)*bcube);
        TVar cosphi = xc[2]/radius;
        TVar xyradius = FADsqrt(xc[0]*xc[0]+xc[1]*xc[1]);
        TVar sinphi = xyradius/radius;
        TVar costheta = xc[0]/xyradius;
        TVar sintheta = xc[1]/xyradius;
        disp[2] = ur*cosphi;
        disp[0] = ur*sinphi*costheta;
        disp[1] = ur*sinphi*sintheta;
    }
    else{
        DebugStop();
    }
}


template<class TVar>
void TElasticity3DAnalytic::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu) const
{
    //    Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    Elast = TVar(fE);
    nu = TVar(fPoisson);
}

//template<>
//void TElasticity3DAnalytic::Elastic(const TPZVec<double> &x, double &Elast, double &nu)
//{
//  Elast = 1000.;
//    Elast = (100. * (1. + 0.3 * sin(10 * M_PI * (x[0] - 0.5)) * cos(10. * M_PI * x[1])));
//    Elast = 1.;
//    nu = 0.;
//}

void TElasticity3DAnalytic::ElasticDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    STATE E = 1.,nu = 0.3;
    DebugStop();
//    Elastic(xstate,E,nu);
    result[0] = E;
    result[1] = nu;
}


template<class TVar>
void TElasticity3DAnalytic::graduxy(const TPZVec<TVar> &x, TPZFMatrix<TVar> &grad) const
{
    int sz = x.size();
#ifdef PZDEBUG
    if(sz != 3)
    {
        DebugStop();
    }
#endif
    TPZManVector<Fad<TVar>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        Fad<TVar> temp = Fad<TVar>(sz,i,x[i]);
        xfad[i] = temp;
        
    }
    TPZManVector<Fad<TVar>,3> result(3);
    uxy(xfad,result);
    grad.Resize(sz,sz);
    for (int i=0; i<sz; i++) {
        for (int j=0; j<sz; j++)
        {
            grad(j,i) = result[i].d(j);
        }
    }
    
}

void TElasticity3DAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    int sz = x.size();
    TPZManVector<Fad<REAL>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        Fad<REAL> temp = Fad<REAL>(sz,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<REAL>,3> result(3);
    uxy(xfad,result);
    gradu.Redim(sz,sz);
    for (int i=0; i<sz; i++) {
        u[i] = result[i].val();
        for (int j=0; j<sz; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

/*
template<>
void TElasticity3DAnalytic::graduxy(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &grad)
{
    TPZManVector<Fad<Fad<REAL> >,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<Fad<REAL> > temp = Fad<Fad<REAL> >(3,Fad<REAL>(3,0.));
        //      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        temp.val()= x[i];
        Fad<REAL> temp3(3,0.);
        for(int j=0; j<3; j++)
        {
            temp.fastAccessDx(j) = temp3;
        }
        Fad<REAL> temp2(3,1.);
        temp.fastAccessDx(i) = temp2;
        //      Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
        //      xfad[i] = temp;
    }
    TPZManVector<Fad<Fad<REAL> >,3> result(2);
    uxy(xfad,result);
    grad.Resize(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            grad(i,j) = result[i].d(j);
        }
    }
}
 */

template<class TVar>
void TElasticity3DAnalytic::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZFNMatrix<9,TVar> grad;
    TVar E, nu;
    Elastic(x, E, nu);
    TVar Fac = E/((TVar)1.+nu)/((TVar(1.)-TVar(2.)*nu));
    graduxy(x,grad);
    //uxy(x,u);
    sigma.Resize(3,3);
    sigma(0,0) = Fac*((TVar(1.)-nu)*grad(0,0)+nu*grad(1,1)+nu*grad(2,2));
    sigma(1,1) = Fac*((TVar(1.)-nu)*grad(1,1)+nu*grad(0,0)+nu*grad(2,2));
    sigma(2,2) = Fac*((TVar(1.)-nu)*grad(2,2)+nu*grad(0,0)+nu*grad(1,1));
    sigma(0,1) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    sigma(0,2) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(0,2)+grad(2,0));
    sigma(2,0) = sigma(0,2);
    sigma(1,2) = E/(TVar(2.)*(TVar(1.)+nu))*(grad(1,2)+grad(2,1));
    sigma(2,1) = sigma(1,2);
}

void TElasticity3DAnalytic::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &tensor) const
{
    TPZManVector<STATE,3> xloc(3,0.);
    for (int i=0; i<3; i++) {
        xloc[i] = x[i];
    }
    Sigma<STATE>(xloc,tensor);
}

/*
template<>
void TElasticity3DAnalytic::Sigma(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma)
{
    TPZFNMatrix<4,Fad<REAL> > grad;
    Fad<REAL>  E, nu;
    Elastic(x, E, nu);
    Fad<REAL>  Fac = E/(Fad<REAL>(1.)+nu)/((Fad<REAL>(1.)-Fad<REAL>(2.)*nu));
    graduxy(x,grad);
    sigma.Resize(2,2);
    sigma(0,0) = Fac*((Fad<REAL>(1.)-nu)*grad(0,0)+nu*grad(1,1));
    sigma(1,1) = Fac*((Fad<REAL>(1.)-nu)*grad(1,1)+nu*grad(0,0));
    sigma(0,1) = E/(Fad<REAL>(2.)*(Fad<REAL>(1.)+nu))*(grad(0,1)+grad(1,0));
    sigma(1,0) = sigma(0,1);
    
}
*/
template
void TElasticity3DAnalytic::Sigma<REAL>(const TPZVec<REAL> &x, TPZFMatrix<REAL> &divsigma) const;

template<class TVar>
void TElasticity3DAnalytic::DivSigma(const TPZVec<REAL> &x, TPZVec<TVar> &divsigma) const
{
    int sz = x.size();
    TPZManVector<Fad<TVar>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        xfad[i] = Fad<TVar>(sz,i,x[i]);
    }
    TPZFNMatrix<9, Fad<TVar> > sigma(3,3);
    Sigma(xfad,sigma);
    for (int i=0; i<3; i++) {
        divsigma[i] = sigma(i,0).dx(0)+sigma(i,1).dx(1)+sigma(i,2).dx(2);
    }
}



template
void TElasticity3DAnalytic::DivSigma<REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &divsigma) const;
template
void TElasticity3DAnalytic::Sigma<Fad<REAL> >(const TPZVec<Fad<REAL> > &x, TPZFMatrix<Fad<REAL> > &sigma) const;


template<class TVar>
void TLaplaceExample1::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    TPZManVector<TVar,3> xloc(x);
    for (int i = 0; i<xloc.size(); i++) {
        xloc[i] -= fCenter[i];
    }
    TVar r2 = 0.;
    for (int i=0; i<fDimension; i++) {
        r2 += xloc[i]*xloc[i];
    }
    TVar r = sqrt(r2);
    disp[0] = xloc[0]*(TVar)(0.);
     
    switch (fExact) {
        case EConst:
            disp[0] += 1.;
            break;
        case EX:
            disp[0] += xloc[0];
            break;
        case ESinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]);
        }
            break;
        case E10SinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]*10.);
        }
            break;
        case ESinDist:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= sin((TVar)M_PI*xloc[i]*1.1);
        }
            break;
        case ECosCos:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= cos((TVar)M_PI*xloc[i]/2.);
        }
            break;
        case EArcTan://(1+0.3sin(10Pi x))*(1+0.5cos(10Pi r)*arctan(100*(r-0.5))
        {
            TVar atanco = (r-(TVar)0.5)*100.;
            TVar freq = 10.;
            TVar mult = (TVar(1)+TVar(0.3)*sin(TVar(M_PI)*xloc[0]*freq))*(TVar(1)+TVar(0.5)*cos(TVar(M_PI)*r*freq));
            disp[0] = atan(atanco)*mult;
        }
            break;
        case EArcTanSingular://5*(0.5*pi + arctang(20(0.25 - r^2))))
        {
            REAL B = 5.;
            if(fDimension==1)
                B *= 0.25;
            else if(fDimension==3)
                B *= 4;
            // Argument value (arc) to compute ArcTangent( arc )
            TVar RCircle = 0.5;
            TVar Force = 20.;
            TVar arc = Force*(RCircle*RCircle-r2);
            TVar Prod = 1.;
            for (int i=0; i<fDimension; i++) {
                Prod *= x[i]*(1.-x[i]);
            }
            TVar temp = 0.5*M_PI + atan(arc);
            disp[0] = B*temp;
        }
            break;
            
            //----
        case ESinMark://(r^(2/3)-r^3)sin(20/3)
        {

            TVar theta = atan2(xloc[1], xloc[0]);//theta=arctan(y/x)

            if (theta < TVar(0.)) theta += 2. * M_PI;

            // Verification to avoid numerical errors when x > 0 and y = 0
            if (xloc[0] > 0 && xloc[1] < 1e-15 && xloc[1] > -1e-15) {
               disp[0] = 0.;
            }
            else {
                TVar factor = pow(r, TVar(2.) / TVar(3.)) - pow(r, TVar(3.));
                disp[0] = factor * ((TVar) (2.) * sin((TVar) (2.) * theta / TVar(3.)));
            }


        }
            break;
            
            //--
     
        case ESinSinDirNonHom: //sin(pi x)sin(pi y)+1/(x+y+1)
        {
            
            disp[0]=(sin((TVar)M_PI*xloc[0]))*sin((TVar)M_PI*xloc[1])+(TVar)(1.)/(xloc[0]+xloc[1]+(TVar)(1.));
            
        }
            
            break;
            
        case ESteklovNonConst://Steklov function for eigenvalue lambda=0.126902 and permeability Omega1=Omega=3, Omega2=Omega4=5
        {
            TVar coefs[] = {1., 0.44721359549995787, 2.3333333333333326,
                -0.7453559924999296, 0.5555555555555556,
                -0.9441175904999111, -0.48148148148148173,
                -2.4017026424997736};
            TVar lambda = 0.53544094560246;
            TVar t = atan2(xloc[1], xloc[0]);
            if(t < TVar(0.)) t += 2.*M_PI;
            
            if((xloc[0] >=TVar(0.)) && (xloc[1] >=TVar(0.))){
               // std::cout<<"1o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                
                disp[0]=pow(r, lambda)*(TVar(coefs[0])*cos(lambda *t) + TVar(coefs[1])*sin(lambda*t) );
               // std::cout<<"valor da funcao no 1o. Q "<<disp[0]<<std::endl;
               // disp[0]=pow(r, lambda)*(cos(lambda *t)+TVar(-1.)*TVar(0.1)*sin(lambda*t));
                
            }
            
            if(( xloc[0] <= TVar(0.)) && (xloc[1] >=TVar(0.))){
               // std::cout<<"2o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[2])*cos(lambda*t) + TVar(coefs[3])*sin(lambda* t));
                
             //    std::cout<<"valor da funcao no 2o. Q "<<disp[0]<<std::endl;
            }
            
            if((xloc[0] <TVar(0.)) && ( xloc[1] <= TVar(0.))){
               // std::cout<<"3o. Q "<<xloc<< " r " << r << " th " << t << std::endl;
                disp[0]= pow(r, lambda)*(TVar(coefs[4] )*cos(lambda*t) + TVar(coefs[5])*sin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(0.882757 )*cos(lambda*t) + TVar(-1.)*TVar(0.480355)*sin(lambda* t));
              //   std::cout<<"valor da funcao no 3o. Q "<<disp[0]<<std::endl;
            }
            if(( xloc[0] >= TVar(0.)) && ( xloc[1] < TVar(0.))){
              //  std::cout<<"4o. Q "<<xloc<< " r " << r << " th " << t << std::endl;

                disp[0]= pow(r, lambda)*(TVar(coefs[6])*cos(lambda*t) +  TVar(coefs[7])*sin(lambda* t));
                
               // std::cout<<"valor da funcao no 4o. Q "<<disp[0]<<std::endl;
                
            }
            
            
        }
            break;
            
            case EGalvisNonConst:
        {
            
            TVar k1 = 2;
            TVar k2 = 5;
            
            
            if((xloc[0] <= TVar(0.)) && (xloc[1] <= TVar(0.))){
                
                TVar u1 = sin(M_PI*(xloc[0]+TVar(1.))/(k1+1));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"3o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if(( xloc[0] > TVar(0.)) && (xloc[1] < TVar(0.))){
             
                
                TVar u1 = sin(M_PI*(k1*xloc[0]+TVar(1.))/(k1 + TVar(1.)));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"4o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if((xloc[0] < TVar(0.)) && ( xloc[1] >= TVar(0.))){
              
                TVar u1 = sin(M_PI*(xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
           //     std::cout<<"2o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            if(( xloc[0] >= TVar(0.)) && ( xloc[1] >= TVar(0.))){
               
                
                TVar u1 = sin(M_PI*(k1*xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"1o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
                
            }
            
            
        }
            break;
            
        case EBoundaryLayer:{
            TVar term1 = xloc[0]*xloc[1]*(1-xloc[0])*(1-xloc[1]);
            TVar term2 = exp(TVar(10.)*xloc[0])*exp(TVar(10.)*xloc[1]);
            TVar factor = TVar(537930);
            
            disp[0] = (term1*term2)/factor;
            
          //  std::cout<<"Pto "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
            
            
        }
            break;
            
        case EBubble:{
            
            disp[0] = xloc[0]*xloc[1]*xloc[2]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1])*(TVar(1.)-xloc[2]);
           // std::cout<<"pto "<<xloc <<" f(x) "<<disp<<std::endl;
        }
            break;
            
        default:
            disp[0] = 0.;
            break;
    }
}

template<>
void TLaplaceExample1::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) const
{
#ifdef PZDEBUG
    
    if(fDimension == -1){
        DebugStop();
    }
    
#endif

    typedef FADFADREAL TVar;
    TPZManVector<TVar,3> xloc(x);
    for (int i = 0; i<xloc.size(); i++) {
        xloc[i] -= fCenter[i];
    }
    TVar r2 = 0.;
    for (int i=0; i<fDimension; i++) {
        r2 += xloc[i]*xloc[i];
    }
    TVar r = FADsqrt(r2);
    disp[0] = (TVar)(0.)*xloc[0];
    switch (fExact) {
        case EConst:
            disp[0] += 1.;
            break;
        case EX:
            disp[0] += xloc[0];
            break;
        case ESinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]);
        }
            break;
        case E10SinSin:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]*10.);
        }
            break;
        case ESinDist:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADsin((TVar)M_PI*xloc[i]*1.1);
        }
            break;
        case ECosCos:
        {
            disp[0] += (TVar)(1.);
            for(int i=0; i<fDimension; i++) disp[0] *= FADcos((TVar)M_PI*xloc[i]/2.);
        }
            break;
        case EArcTan:
        {
            TVar atanco = (r-(TVar)0.5)*100.;
            TVar freq = 10.;
            TVar mult = (TVar(1)+TVar(0.3)*FADsin(TVar(M_PI)*xloc[0]*freq))*(TVar(1)+TVar(0.5)*FADcos(TVar(M_PI)*r*freq));
            disp[0] = FADatan(atanco)*mult;
        }
            break;
        case EArcTanSingular:
        {
            REAL B = 5.;
            if(fDimension==1)
                B *= 0.25;
            else if(fDimension==3)
                B *= 4;
            // Argument value (arc) to compute ArcTangent( arc )
            TVar RCircle = 0.5;
            TVar Force = 20.;
            TVar arc = Force*(RCircle*RCircle-r2);
            TVar Prod = 1.;
            for (int i=0; i<fDimension; i++) {
                Prod *= x[i]*(1.-x[i]);
            }
            TVar temp = 0.5*M_PI + FADatan(arc);
            disp[0] = B*temp;
        }
            break;
        case ESinSinDirNonHom: //sin(pi x)sin(pi y)+1/(x+y+1)
        {
            
            disp[0]=(FADsin((TVar)M_PI*xloc[0]))*FADsin((TVar)M_PI*xloc[1])+(TVar)(1.)/(xloc[0]+xloc[1]+(TVar)(1.));
            
        }
            
            break;
        case ESinMark://(r^(2/3)-r^3)sin(20/3)
        {
            
            TVar theta=FADatan2(xloc[1],xloc[0]);//theta=atan(y/x)
            if( theta < TVar(0.)) theta += 2.*M_PI;
            
            TVar factor=pow(r,TVar (2.)/TVar (3.))-pow(r,TVar (3.));
            disp[0]= factor*((TVar)(2.)*FADsin((TVar)(2.)*theta/TVar(3.)));
            
        }
            break;
        case ESteklovNonConst://Steklov function for eigenvalue lambda=0.126902 and permeability Omega1=Omega=3=100, Omega2=Omega4=1
        {
            
            TVar coefs[] = {1., 0.44721359549995787, 2.3333333333333326,
                -0.7453559924999296, 0.5555555555555556,
                -0.9441175904999111, -0.48148148148148173,
                -2.4017026424997736};
            TVar lambda = 0.53544094560246;
            TVar t = FADatan2(xloc[1], xloc[0]);
            if(t < TVar(0.)) t += TVar(2.*M_PI);

            if((xloc[0] >= TVar(0.)) && (xloc[1] >= TVar(0.))){
              //  std::cout<<"1o. Q"<<xloc<<std::endl;
                
                disp[0]=pow(r, lambda)*(TVar(coefs[0])*FADcos(lambda *t) + TVar(coefs[1])*FADsin(lambda*t));
              //  std::cout<<"valor da funcao no 1o. Q "<<disp[0]<<std::endl;
                // disp[0]=pow(r, lambda)*(cos(lambda *t)+TVar(-1.)*TVar(0.1)*sin(lambda*t));
                
            }
            
            if(( xloc[0] < TVar(0)) && (xloc[1] >TVar(0.))){
               // std::cout<<"2o. Q"<<xloc<<std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[2])*FADcos(lambda*t) + TVar(coefs[3])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(2.9604)*cos(lambda*t) +TVar(-1.)* TVar(9.60396)*sin(lambda* t));
                //std::cout<<"valor da funcao no 2o. Q "<<disp[0]<<std::endl;
            }
            
            if((xloc[0] < TVar(0.)) && ( xloc[1] <= TVar(0.))){
             //   std::cout<<"3o. Q"<<xloc<<std::endl;
                disp[0]= pow(r, lambda)*(TVar(coefs[4] )*FADcos(lambda*t) + TVar(coefs[5])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(0.882757 )*cos(lambda*t) + TVar(-1.)*TVar(0.480355)*sin(lambda* t));
               // std::cout<<"valor da funcao no 3o. Q "<<disp[0]<<std::endl;
            }
            if(( xloc[0] >= TVar(0.)) && ( xloc[1] < TVar(0.))){
               // std::cout<<"4o. Q"<<xloc<<std::endl;
                
                disp[0]= pow(r, lambda)*(TVar(coefs[6])*FADcos(lambda*t) +  TVar(coefs[7])*FADsin(lambda* t));
                //disp[0]= pow(r, lambda)*(TVar(-1.)*TVar(6.45646)*cos(lambda*t) +  TVar(7.70156 )*sin(lambda* t));
              //  std::cout<<"valor da funcao no 4o. Q "<<disp[0]<<std::endl;
                
            }

            
            
        }
            break;
            
        case EGalvisNonConst:
        {
            
            TVar k1 = 2;
            TVar k2 = 5;
            
            
            if((xloc[0] <= TVar(0.)) && (xloc[1] <= TVar(0.))){
                
                TVar u1 = FADsin(M_PI*(xloc[0]+TVar(1.))/(k1+1));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"3o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if(( xloc[0] > TVar(0.)) && (xloc[1] < TVar(0.))){
                
                
                TVar u1 = FADsin(M_PI*(k1*xloc[0]+TVar(1.))/(k1 + TVar(1.)));
                TVar u2 = TVar(4.)*(xloc[1]+TVar(1.))*(k2-xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
            //    std::cout<<"4o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            
            if((xloc[0] < TVar(0.)) && ( xloc[1] >= TVar(0.))){
                
                TVar u1 = FADsin(M_PI*(xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
              //  std::cout<<"2o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
            }
            if(( xloc[0] >= TVar(0.)) && ( xloc[1] >= TVar(0.))){
                
                
                TVar u1 = FADsin(M_PI*(k1*xloc[0]+TVar(1.))/(k1+TVar(1.)));
                TVar u2 = TVar(4.)*(k2*xloc[1]+TVar(1.))*(k2-k2*xloc[1])/((k2+TVar(1.))*(k2+TVar(1.)));
                
                disp[0]=u1*u2;
                //std::cout<<"1o. Q "<<xloc<<" valor da funcao "<<disp[0]<<std::endl;
                
                
            }
        }
            break;
            
        case EBoundaryLayer:{
            TVar term1 = xloc[0]*xloc[1]*(1-xloc[0])*(1-xloc[1]);
            TVar term2 = FADexp(TVar(10.)*xloc[0])*FADexp(TVar(10.)*xloc[1]);
            TVar factor = TVar(537930);
            
            disp[0] = (term1*term2)/factor;

            
        }
            break;
            
        case EBubble:{
            
            disp[0] = xloc[0]*xloc[1]*xloc[2]*(TVar(1.)-xloc[0])*(TVar(1.)-xloc[1])*(TVar(1.)-xloc[2]);
         //   std::cout<<"pto "<<xloc <<" f(x) "<<disp<<std::endl;
            
        }
            break;
            
        default:
            disp[0] = xloc[0]*0.;
            break;
    }

}

template<class TVar>
void TLaplaceExample1::Permeability(const TPZVec<TVar> &x, TVar &Perm)
{
    Perm = (TVar)(1.);
}

void TLaplaceExample1::PermeabilityDummy(const TPZVec<REAL> &x, TPZVec<STATE> &result, TPZFMatrix<STATE> &deriv)
{
    TPZManVector<STATE,3> xloc(x.size());
    for (unsigned int i = 0; i<xloc.size(); ++i) {
        xloc[i] = x[i];
    }
    STATE Perm;
    Permeability(xloc, Perm);
    deriv.Zero();
    deriv(0,0) = Perm;
    deriv(1,1) = Perm;
    deriv(2,2) = Perm;
    deriv(3,0) = 1./Perm;
    deriv(4,1) = 1./Perm;
    deriv(5,2) = 1./Perm;
}

template<class TVar>
void TLaplaceExample1::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<TVar> temp = Fad<TVar>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(3);
    for (int i=0; i<3; i++)
    {
            grad[i] = result[0].d(i);
    }
}

void TLaplaceExample1::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<REAL> temp = Fad<REAL>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<REAL>,3> result(1);
    uxy(xfad,result);
    gradu.Redim(3,1);
    
    
    u[0] = result[0].val();
    for (int i=0; i<3; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }

    
}

template<class TVar>
void TLaplaceExample1::SigmaLoc(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZManVector<TVar,3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x,grad);
    sigma.Resize(3,1);
    sigma(0) = -Perm*grad[0];
    sigma(1) = -Perm*grad[1];
    sigma(2) = -Perm*grad[2];

}

template<class TVar>
void TLaplaceExample1::DivSigma(const TPZVec<TVar> &x, TVar &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        xfad[i] = Fad<TVar>(3,i,x[i]);
    }
    TPZFNMatrix<3,Fad<TVar> > sigma(3,1);
    SigmaLoc(xfad,sigma);
    divsigma = sigma(0).dx(0)+sigma(1).dx(1)+sigma(2).dx(2);
    
}

template
void TLaplaceExample1::SigmaLoc(const TPZVec<STATE> &x, TPZFMatrix<STATE> &sigma) const;

template
void TLaplaceExample1::DivSigma<REAL>(const TPZVec<REAL> &x, REAL &divsigma) const;



//ExactFunc *Exact();


template<class TVar>
void TLaplaceExampleTimeDependent::uxy(const TPZVec<TVar> &x, TPZVec<TVar> &disp) const
{
    switch(fProblemType)
    {
        case ELinear:
            disp[0] = x[0];
            break;
        case ESin:
            disp[0] = sin(M_PI*x[0])*sin(M_PI*x[1])*exp(-2.*M_PI*M_PI*fTime);
            break;
        case ECos:
            disp[0] = cos(M_PI_2*x[0])*cos(M_PI_2*x[1])*exp(-M_PI*M_PI_2*fTime);
            break;
        default:
            DebugStop();
    }
}

template<>
void TLaplaceExampleTimeDependent::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &disp) const
{
    switch(fProblemType)
    {
        case ELinear:
            disp[0] = x[0];
            break;
        case ESin:
            disp[0] = FADsin(M_PI*x[0])*FADsin(M_PI*x[1])*FADexp(-2.*M_PI*M_PI*fTime);
            break;
        case ECos:
            disp[0] = FADcos(M_PI_2*x[0])*FADcos(M_PI_2*x[1])*FADexp(-M_PI*M_PI_2*fTime);
            break;
        default:
            DebugStop();
    }

}

template<class TVar>
void TLaplaceExampleTimeDependent::Permeability(const TPZVec<TVar> &x, TVar &Perm) const
{
    Perm = (TVar)(fK);
}


template<class TVar>
void TLaplaceExampleTimeDependent::graduxy(const TPZVec<TVar> &x, TPZVec<TVar> &grad) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<TVar> temp = Fad<TVar>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<TVar>,3> result(1);
    uxy(xfad,result);
    grad.resize(2);
    for (int i=0; i<2; i++)
    {
        grad[i] = result[0].d(i);
    }
}

void TLaplaceExampleTimeDependent::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,x[i]);
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,1);
    u[0] = result[0].val();
    for (int i=0; i<2; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

template<class TVar>
void TLaplaceExampleTimeDependent::Sigma(const TPZVec<TVar> &x, TPZFMatrix<TVar> &sigma) const
{
    TPZManVector<TVar,3> grad;
    TVar Perm;
    Permeability(x, Perm);
    graduxy(x,grad);
    sigma.Resize(2,1);
    sigma(0) = -Perm*grad[0];
    sigma(1) = -Perm*grad[1];
    
}


void TLaplaceExampleTimeDependent::Sigma(const TPZVec<REAL> &x, TPZFMatrix<STATE> &sigma) const
{
    typedef STATE TVar;
    TPZManVector<TVar,3> grad, xst(3);
    for(int i=0; i<3; i++) xst[i] = x[i];
    TVar Perm;
    Permeability(xst, Perm);
    graduxy(xst,grad);
    sigma.Resize(2,1);
    sigma(0) = -Perm*grad[0];
    sigma(1) = -Perm*grad[1];
    
}

template<class TVar>
void TLaplaceExampleTimeDependent::DivSigma(const TPZVec<TVar> &x, TVar &divsigma) const
{
    TPZManVector<Fad<TVar>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        xfad[i] = Fad<TVar>(2,i,x[i]);
    }
    TPZFNMatrix<3, Fad<TVar> > sigma(2,1);
    Sigma(xfad,sigma);
    divsigma = sigma(0).dx(0)+sigma(1).dx(1);
    
}


template
void TLaplaceExampleTimeDependent::DivSigma(const TPZVec<REAL> &x, REAL &divsigma) const;


template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::uxy(const TPZVec<TVar1> &x, TPZVec<TVar2> &flux) const
{
    TVar1 x1 = x[0];
    TVar1 x2 = x[1];
    
    switch(fProblemType)
    {
        case EStokes0:
            flux[0] = -0.1*x2*x2+0.2*x2;
            flux[1] = 0.;
            break;
        case EStokesLimit1:
            flux[0] = -1.*sin(x1)*sin(x2);;
            flux[1] = -1.*cos(x1)*cos(x2);
            break;
        case EBrinkman1:
            flux[0] = 0.;
              break;
        case EDarcyLimit1:
            flux[0] = 0.;
            break;
        default:
            DebugStop();
    }
}


template<>
void TStokes2DAnalytic::uxy(const TPZVec<FADFADREAL > &x, TPZVec<FADFADREAL > &flux) const
{
    FADFADREAL x1 = x[0];
    FADFADREAL x2 = x[1];

    switch(fProblemType)
    {
        case EStokes0:
            flux[0] = -0.1*x2*x2+0.2*x2;
            flux[1] = 0.;
            break;
        case EStokesLimit1:
            flux[0] = -1.*FADsin(x1)*FADsin(x2);;
            flux[1] = -1.*FADcos(x1)*FADcos(x2);
            break;
        case EBrinkman1:
            flux[0] = 0.;
            break;
        case EDarcyLimit1:
            flux[0] = 0.;
            break;
        default:
            DebugStop();
    }
    
}

template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::pressure(const TPZVec<TVar1> &x, TVar2 &p) const
{
    TVar1 x1 = x[0];
    TVar1 x2 = x[1];
    
    switch(fProblemType)
    {
        case EStokes0:
            p = 1.-0.2*x1;
            break;
        case EStokesLimit1:
            p = cos(x1)*sin(x2);
            break;
        case EBrinkman1:
            p = 0.;
            break;
        case EDarcyLimit1:
            p = 0.;
            break;
        default:
            DebugStop();
    }
}

template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::graduxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &gradu) const
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<REAL> temp = Fad<REAL>(2,i,shapeFAD::val(x[i]));
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<REAL>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
}

template<>
void TStokes2DAnalytic::graduxy(const TPZVec<std::complex<double> > &x, TPZFMatrix<std::complex<double> > &gradu) const
{
    TPZManVector<Fad<std::complex<double> >,3> xfad(x.size());
    for(int i=0; i<2; i++)
    {
        Fad<std::complex<double>> temp = Fad<std::complex<double>>(2,i,shapeFAD::val(x[i]));
        xfad[i] = temp;
    }
    xfad[2] = x[2];
    TPZManVector<Fad<std::complex<double>>,3> result(2);
    uxy(xfad,result);
    gradu.Redim(2,2);
    for (int i=0; i<2; i++) {
        for (int j=0; j<2; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
}

template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::Duxy(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &Du) const
{
    TPZFMatrix<TVar2> grad(3,3,0.), gradT(3,3,0.);
    graduxy(x,grad);
    grad.Transpose(&gradT);
    Du = (grad+gradT)*(TVar2)0.5;
}

template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::Sigma(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &sigma) const
{
    SigmaLoc(x, sigma);
}


template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::SigmaLoc(const TPZVec<TVar1> &x, TPZFMatrix<TVar2> &sigma) const
{
    TPZFMatrix<TVar2> Du, pIdentity(sigma.Rows(),sigma.Cols());
    TVar2 p=0.;
    Duxy(x,Du);
    for(int i=0; i<Du.Rows(); i++)
    {
        for(int j=0; j<Du.Cols(); j++)
        {
            Du(i,j) *= 2.*fvisco;
        }
    }
    pressure(x, p);
    for (int i=0; i< pIdentity.Rows(); i++) {
        pIdentity(i,i) = p;
    }
    sigma = Du-pIdentity;
}

template<typename TVar1, typename TVar2>
void TStokes2DAnalytic::DivSigma(const TPZVec<TVar1> &x, TPZVec<TVar2> &divsigma) const
{
    int sz = x.size();
    TPZManVector<Fad<TVar2>,3> xfad(sz);
    for(int i=0; i<sz; i++)
    {
        xfad[i] = Fad<TVar2>(sz,i,x[i]);
    }
    TPZFNMatrix<9, Fad<TVar2> > sigma(3,3);
    Sigma(xfad,sigma);
    for (int i=0; i<3; i++) {
        divsigma[i] = sigma(i,0).dx(0)+sigma(i,1).dx(1)+sigma(i,2).dx(2);
    }
}

void TStokes2DAnalytic::Solution(const TPZVec<REAL> &x, TPZVec<STATE> &u, TPZFMatrix<STATE> &gradu) const
{
    TPZManVector<Fad<REAL>,3> xfad(x.size());
    for(int i=0; i<3; i++)
    {
        Fad<REAL> temp = Fad<REAL>(3,i,x[i]);
        xfad[i] = temp;
    }
    TPZManVector<Fad<REAL>,3> result(1);
    uxy(xfad,result);
    gradu.Redim(3,1);
    u[0] = result[0].val();
    for (int i=0; i<3; i++) {
        for (int j=0; j<1; j++)
        {
            gradu(i,j) = result[j].d(i);
        }
    }
    
}

#endif
