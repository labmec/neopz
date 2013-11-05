 #include "poroelastoplastic.h"
//#include "pzskylstrmatrix.h"
//#include "TPZReadGIDGrid.h"
#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
//#include "pzcompelwithmem.h"
//#include "pzelastoplastic.h"
//#include "pzporous.h"
#include "TPZLadeKim.h"
//#include "TPZSandlerDimaggio.h"
//#include "TPZYCDruckerPrager.h"
//#include "TPZThermoForceA.h"
//#include "TPZElasticResponse.h"
//#include "pzelastoplasticanalysis.h"
#include "pzmat2dlin.h"
#include "pzporoanalysis.h"
//#include "TPZTensor.h"
//#include "BrazilianTestGeoMesh.h"
//#include "pzelast3d.h"
//#include "pzcompelpostproc.h"
//#include "pzpostprocmat.h"
//#include "pzpostprocanalysis.h"
//#include "pzblockdiag.h"
//#include "TPZSpStructMatrix.h"
//#include "pzbdstrmatrix.h"
//#include "pzstepsolver.h"
#include "pzbfilestream.h"
#include <sstream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>
//#include "pzelastoplastic.h"
//#include "pzporous.h"
//#include "TPZThermoForceA.h"
//#include "TPZElasticResponse.h"
//#include "pzelastoplasticanalysis.h"
//#include "pzanalysis.h"
//#include "pzskylstrmatrix.h"
//#include "TPZTensor.h"
//#include "pzcompelpostproc.h"
//#include "pzpostprocmat.h"
//#include "pzpostprocanalysis.h"
//#include "TPZYCVonMises.h"
//#include "TPZVonMises.h"
//#include "pzfstrmatrix.h"
//#include "pzbndmat.h"
//#include "pzgeoquad.h"
//#include "TPZGeoCube.h"
//#include "pzgeotetrahedra.h"
//#include "pzgeopyramid.h"
//#include "tpzgeoelrefpattern.h"
//#include "pzbndcond.h"
//#include "pzstepsolver.h"
//#include "TPZTensor.h"
//#include "TPZYCMohrCoulomb.h"
//#include "TPZMohrCoulomb.h"
//#include "TPZDruckerPrager.h"
#include "pzelasmat.h"
//#include "pzelastoplastic2D.h"
//#include "tpzycvonmisescombtresca.h"
//#include "TPZMohrCoulombNeto.h"
//#include "TPZSandlerDimaggio.h"
#include "TPZVTKGeoMesh.h"
#include "BrazilianTestGeoMesh.h"

#include "pzsandlerextPV.h"

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

using namespace pzshape; // needed for TPZShapeCube and related classes

#include <math.h>

#define MACOS
#ifdef MACOS

#include <iostream>
#include <math.h>
#include <signal.h>
#include <fenv.h>
#include <xmmintrin.h>

#define ENABLE_FPO_EXCEPTIONS _MM_SET_EXCEPTION_MASK(_MM_GET_EXCEPTION_MASK() & ~_MM_MASK_INVALID);


#define DECLARE_FPO_HANDLER_FUNC void InvalidFPOHandler(int signo) {\
switch(signo) {\
case SIGFPE: std::cout << "ERROR : Invalid Arithmetic operation." << std::endl; break;\
}\
exit(signo);\
}

#define ATTACH_FPO_SIGNAL struct sigaction act = {};\
act.sa_handler = InvalidFPOHandler;\
sigaction(SIGFPE, &act, NULL);


DECLARE_FPO_HANDLER_FUNC;
#endif

#include "WellBoreAnalysis.h"
#include "pzsandlerextPV.h"

int main()
{
    
//    // D ->  0.67, R -> 2.5(*5.3307429045*), A ->  0.25, B ->  0.67, C -> \
//    // 0.18, W ->  0.066, \[Psi] -> 1, NN -> 0, \[Phi] -> 0, G -> 40, K -> \
//    // 66.6667, d -> 0.67, L0 -> -0.123148527544659
//    //    REAL A, REAL B,REAL C, REAL D,REAL K,REAL G,REAL W,REAL R,REAL Phi,REAL N,REAL Psi
//    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
//    REAL val = materialmodel.F(0.1,0);
//    REAL val2 = materialmodel.X(0.13);
//    REAL val3 = materialmodel.EpsEqX((val2));
//    REAL val4 = materialmodel.EpsEqk(0.13);
//    TPZTensor<REAL> a,aa;
//    a.XX()=0.1;    a.XY()=0;    a.XZ()=0;    a.YY()=-0.56;    a.YZ()=0;    a.ZZ()=-0.3;
//    a.XX()=0.12;    a.XY()=0;    a.XZ()=0;    a.YY()=0.18;    a.YZ()=0;    a.ZZ()=0.05;
//    TPZTensor<REAL>::TPZDecomposed Tensor(a),proj;
//    REAL k0,epsp=0.;
//    materialmodel.Firstk(epsp,k0);
////    REAL res = materialmodel.ResL(Tensor, M_PI/3.,M_PI, 0.12, k0);
//    TPZManVector<REAL> solf1,solf2;
//    solf1=materialmodel.F1Cyl(-0.1,M_PI);
//    solf2 = materialmodel.F2Cyl(M_PI/3,-0.11,k0);
////    
////    REAL distf1 = materialmodel.DistF1(Tensor,-0.1,M_PI);
////    REAL distf2 = materialmodel.DistF2(Tensor,M_PI/3,M_PI,k0);
////    TPZFMatrix<REAL> d2distf21=materialmodel.DDistFunc2(Tensor,M_PI/3,M_PI,0.12,k0);
////    cout << "\n DDistFunc2 = "<<d2distf21<<endl;
////    TPZFMatrix<REAL> ddistf1=materialmodel.DDistFunc1(Tensor,-0.1,M_PI);
////    d2distf21=materialmodel.D2DistFunc2(Tensor,M_PI/3,M_PI,0.133045);
////    cout << "\n D2DistFunc2 = "<<d2distf21<<endl;
////    cout << "\n resf1 = "<<ddistf1<<endl;
//    TPZVec<REAL> vec;
//    aa.XX()=0.12;    aa.XY()=0;    aa.XZ()=0;    aa.YY()=0.18;    aa.YZ()=0;    aa.ZZ()=0.05;
//    TPZTensor<REAL>::TPZDecomposed Tensor2(aa);
//    materialmodel.Firstk(epsp,k0);
//    TPZVec<REAL> yield(2);
//    materialmodel.ProjectF2(Tensor, proj,k0);
//    solf1=materialmodel.FromPrincipalToHWCyl(proj);
//    cout << "\n solf1 = "<<solf1<<endl;
//    REAL beta=solf1[2];
//    
//    materialmodel.ProjectF1(Tensor, proj);
//    solf1=materialmodel.FromPrincipalToHWCyl(proj);
//    beta=solf1[2];
//    materialmodel.YieldFunction(proj, yield, k0);
//    cout << "\n yield = "<<yield<<endl;
//    
//    
//    materialmodel.Firstk(epsp,k0);
//    materialmodel.ProjectRing(Tensor2, proj,k0);
//    solf1=materialmodel.FromPrincipalToHWCyl(proj);
//    cout << "\n solf1 = "<<solf1<<endl;
//    beta=solf1[2];
//    materialmodel.YieldFunction(proj, yield, k0);
//    cout << "\n yield = "<<yield<<endl;
//    //materialmodel.ProjectRing(proj,proj,k0);
//    //solf1=materialmodel.FromPrincipalToHWCyl(proj);
//    //cout << "\n solf1 = "<<solf1<<endl;
    
/*    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
    TPZPlasticState<REAL> plasticstate;
    TPZManVector<STATE,3> epst(3);
//    epst.XX()=-0.001;
//    epst.XY()=0.;
//    epst.XZ()=0.;
//    epst.YY()=-0.0056;
//    epst.YZ()=0.;
//    epst.ZZ()=-0.003;
    epst[0]=-0.001;
    epst[1]=0.;
    epst[2]=0.;
    

   // materialmodel.ApplyStrainComputeSigma(plasticstate,stress);
    
    ofstream outfiletxty("FIGURA12amodeloDiogo.txt");
    TPZManVector<STATE,3> deltaeps(3),eps(3),sigma(3),deltasigma(3);
    

    deltaeps[0]= -0.00135;
    deltaeps[1]=0;
    deltaeps[2]=0;
    eps=deltaeps;
    STATE kprev=0.,epspv1=0.,sig1,eps1;
    materialmodel.Firstk(epspv1,kprev);
    TPZManVector<STATE,3> Tensor2(3);
    Tensor2[0]=-0.12;
    Tensor2[1]=-0.18;
    Tensor2[2]=-0.25;
    materialmodel.CheckCoordinateTransformation(Tensor2);
    
    TPZManVector<STATE,2> yield(2);
    materialmodel.YieldFunction(Tensor2,kprev,yield);
    cout << "yield = "<<yield <<endl;
    
    // perform convergence verifications
    STATE theta = M_PI/5.;
    STATE xi = (Tensor2[0]+Tensor2[1]+Tensor2[2])/sqrt(3.);
    STATE beta = M_PI/3.;
    STATE k = -0.3;
    TPZManVector<STATE> xnorm,errnorm,converge;
    
    materialmodel.TaylorCheckDistF1(Tensor2, xi, beta, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    
    materialmodel.TaylorCheckDDistF1DSigtrial(Tensor2, xi, beta, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    
    materialmodel.TaylorCheckDistF2(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckDDistF2(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckDDistF2DSigtrial(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "xnorm " << xnorm << endl;
    cout << "errnorm " << errnorm << endl;
    cout << "convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckDF1Cart(xi, beta, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckDF2Cart(theta, beta, k, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckProjectSigma(Tensor2, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    
    materialmodel.TaylorCheckParamF1Sigtrial(Tensor2,kprev,xnorm,errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    
    materialmodel.TaylorCheckProjectF1(Tensor2, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    
    materialmodel.TaylorCheckProjectF2(Tensor2, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
  */  
/*
    
   // ApplyStrainComputeSigma(TPZVec<STATE> &epst,TPZVec<STATE> &epsp,STATE & kprev,TPZVec<STATE> &epspnext,TPZVec<STATE> &stressnext,STATE & knext) const
    TPZManVector<STATE,3> epsPnext(3),stresnext(3),epsT(3),epsP(3),deleps(3),STRE(3),STRENEXT(3),epssol(3);
    STATE knext,kpr;
    
    ofstream outfile("FIGURA_12.txt");
    
    
    deleps[0]= -0.00135;
    deleps[1]=0;
    deleps[2]=0;
    epsT[0]=-0.0000001;
    epsT[1]=0.;
    epsT[2]=0.;
    
    epsP[0]=0.;
    epsP[1]=0.;
    epsP[2]=0.;
    
    
    //STATE epv=0.;
    //materialmodel.Firstk(epv,kpr);
    kpr=0.13;
    
    for(int i=0;i<60;i++)
    {
        materialmodel.ApplyStrainComputeSigma(epsT,epsP,kpr,epsPnext,stresnext,knext);
        if(i==50)
        {
            deleps[0]*=-1;
        }
        
        sig1 = stresnext[0];
        eps1 =  epsT[0];
        outfile << -eps1<< " " << -sig1 << "\n";
        
        
        epsT[0]+=deleps[0];
        epsT[1]+=deleps[1];
        epsT[2]+=deleps[2];
        
        epsP[0]=epsPnext[0];
        epsP[1]=epsPnext[1];
        epsP[2]=epsPnext[2];
        
        kpr=knext;
    }
    
    ofstream outfile2("FIGURA_11a.txt");
    ofstream outfile3("FF.txt");
    
    
    deleps[0]= -0.00135;
    deleps[1]=0;
    deleps[2]=0;
    epsT[0]=-0.0000001;
    epsT[1]=0.;
    epsT[2]=0.;
    
    epsP[0]=0.;
    epsP[1]=0.;
    epsP[2]=0.;
    
    STATE I1,sqrtj2,I1tr;
    //materialmodel.Firstk(epv,kpr);
    kpr=0.13;
    
    for(int i=0;i<60;i++)
    {
        materialmodel.ApplyStrainComputeSigma(epsT,epsP,kpr,epsPnext,stresnext,knext);
        if(i==50)
        {
            deleps[0]*=-1;
        }
        
        materialmodel.ComputeI1(stresnext,I1);
        materialmodel.ComputeJ2(stresnext,sqrtj2);
        sqrtj2=sqrt(sqrtj2);
        outfile2 << -I1 << " " << sqrtj2 << "\n";
      

        materialmodel.ApplyStrainComputeElasticStress(STRE, epsT);
        
        materialmodel.ComputeI1(STRE,I1tr);
        STATE F = materialmodel.F(I1tr,0);
        outfile3 << -I1tr << " " << F << "\n";
        
        
        epsT[0]+=deleps[0];
        epsT[1]+=deleps[1];
        epsT[2]+=deleps[2];
        
        epsP[0]=epsPnext[0];
        epsP[1]=epsPnext[1];
        epsP[2]=epsPnext[2];
        
        kpr=knext;
    }
    
*/
    
    
    //void TPZSandlerExtended::ProjectSigmaDep(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &sigproj,STATE &kproj, TPZFMatrix<STATE> &GradSigma) const

    // ApplyStrainComputeSigma(TPZVec<STATE> &epst,TPZVec<STATE> &epsp,STATE & kprev,TPZVec<STATE> &epspnext,TPZVec<STATE> &stressnext,STATE & knext) const
    
    
    
    TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3);
    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
    
    ofstream outfile("FIGURA_12x.txt");
    
    deleps[0]= -0.00135;
    
    for (int k=0; k<3;k++) {
        deltaepsP[k]=0.;
        epsT[0]=0.;
        epsPnext[0]=0.;
        epssol[0]=0.;
    }
    
    STATE kproj,kprev,epv=0.;
    TPZFMatrix<STATE> GradSigma;
    
    kproj=0.;
    kprev=0.13;
    materialmodel.Firstk(epv,kprev);
    for(int i=0;i<65;i++)
    {

        for (int k=0; k<3;k++) {
            epssol[k]=epsT[k]-epsPnext[k];
        }
        
        materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial);
        materialmodel.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
        outfile << -epsT[0]<< " " << -sigproj[0] << "\n";

        if(i==50)
        {
            deleps[0]*=-1;
        }
        
        for (int k=0; k<3;k++) {
            deltasigma[k]=sigtrial[k]-sigproj[k];
        }
        
        materialmodel.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);

        for (int k=0; k<3;k++) {
            epsPnext[k]+=deltaepsP[k];
            epsT[k]+=deleps[k];
        }
        kprev=kproj;
        
    }
    
    
    return 0;
}

