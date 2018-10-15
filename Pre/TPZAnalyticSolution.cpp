
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
            disp[0] = x[0]*(1.-fPoisson*fPoisson)/fE;
            disp[1] = -x[1]*(1.+fPoisson)*fPoisson/fE;
        }
        else
        {
            disp[0] = x[0]/fE;
            disp[1] = -x[1]*fPoisson/fE;
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
    else if(fProblemType==EBend)
    {
        typedef TVar FADFADREAL;
        TVar poiss = fPoisson;
        TVar elast = fE;
        if(fPlaneStress == 0)
        {
            poiss = poiss/(1.-poiss);
            elast /= (1-fPoisson*fPoisson);
        }
        disp[0] = 5.*x[0]*x[1]/elast;
        disp[1] = (-poiss*5.*x[1]*x[1]/2.-5.*x[0]*x[0]/2.)/elast;
        
    }
    else if(fProblemType == ELoadedBeam)
    {
        TVar Est,nust,G;
        REAL MI = 5, h = 1.;
        G = fE/(2.*(1.+fPoisson));
        if (fPlaneStress == 0) {
            Est = fE/((1.-fPoisson*fPoisson));
            nust = fPoisson/(1-fPoisson);
//            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
//            nust = fPoisson/(1.+fPoisson);
        }
        else
        {
            Est = fE;
            nust = fPoisson;
        }
        disp[0] = MI*h*h*x[1]/(2.*G)+ MI*x[0]*x[0]*x[1]/(2.*Est)-MI *x[1]*x[1]*x[1]/(6.*G)+MI*nust*x[1]*x[1]*x[1]/(6.*Est);
        disp[1] = -MI*x[0]*x[0]*x[0]/(6*Est)-MI*nust*x[0]*x[1]*x[1]/(2.*Est);
    }
    else if(fProblemType == ESquareRoot)
    {
        TVar Est,nust,G, kappa;
        TVar theta = FADatan2(x[1],x[0]);
        TVar r = FADsqrt(x[0]*x[0]+x[1]*x[1]);
        G = fE/(2.*(1.+fPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
            //            nust = fPoisson/(1.+fPoisson);
            Est = fE/((1.-fPoisson*fPoisson));
            nust = fPoisson/(1-fPoisson);
            kappa = 3.-4.*fPoisson;
        }
        else
        {
            Est = fE;
            nust = fPoisson;
            kappa = (3.-fPoisson)/(1+fPoisson);
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
        G = fE/(2.*(1.+fPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
            //            nust = fPoisson/(1.+fPoisson);
            Est = fE/((1.-fPoisson*fPoisson));
            nust = fPoisson/(1-fPoisson);
            kappa = 3.-4.*fPoisson;
        }
        else
        {
            Est = fE;
            nust = fPoisson;
            kappa = (3.-fPoisson)/(1+fPoisson);
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
        G = fE/(2.*(1.+fPoisson));
        if (fPlaneStress == 0) {
            //            Est = (1.+2.*fPoisson)/((1+fPoisson)*(1.+fPoisson))*fE;
            //            nust = fPoisson/(1.+fPoisson);
            Est = fE/((1.-fPoisson*fPoisson));
            nust = fPoisson/(1-fPoisson);
            kappa = 3.-4.*fPoisson;
        }
        else
        {
            Est = fE;
            nust = fPoisson;
            kappa = (3.-fPoisson)/(1+fPoisson);
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

template
void TElasticity2DAnalytic::uxy<REAL,REAL>(const TPZVec<REAL> &x, TPZVec<REAL> &disp) const;

#if (REAL != STATE)

template
void TElasticity2DAnalytic::uxy<REAL,STATE>(const TPZVec<REAL> &x, TPZVec<STATE> &disp) const;

template
void TElasticity2DAnalytic::uxy<STATE,STATE>(const TPZVec<STATE> &x, TPZVec<STATE> &disp) const;

#endif

template<class TVar>
void TElasticity2DAnalytic::Elastic(const TPZVec<TVar> &x, TVar &Elast, TVar &nu) const
{
//    Elast = (TVar(100.) * (TVar(1.) + TVar(0.3) * sin(TVar(10 * M_PI) * (x[0] - TVar(0.5))) * cos(TVar(10. * M_PI) * x[1])));
    Elast = fE;
    nu = TVar(fPoisson);
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
    TPZManVector<STATE> xstate(x.size());
    for (int i=0; i<xstate.size(); i++) {
        xstate[i] = x[i];
    }
    STATE E = 1,nu = 0.3;
//    Elastic(xstate,E,nu);
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
void TElasticity2DAnalytic::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma) const
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
void TElasticity3DAnalytic::DivSigma(const TPZVec<TVar> &x, TPZVec<TVar> &divsigma) const
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
        case EArcTan:
        {
            TVar atanco = (r-(TVar)0.5)*100.;
            TVar freq = 10.;
            TVar mult = (TVar(1)+TVar(0.3)*sin(TVar(M_PI)*xloc[0]*freq))*(TVar(1)+TVar(0.5)*cos(TVar(M_PI)*r*freq));
            disp[0] = atan(atanco)*mult;
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
            TVar temp = 0.5*M_PI + atan(arc);
            disp[0] = B*temp;
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

#endif
