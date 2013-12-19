#include "pzmat2dlin.h"
#include "pzporoanalysis.h"
#include "pzbfilestream.h"
#include <sstream>
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif
#include <iostream>
#include <cstdlib>
#include "pzelasmat.h"
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

#include "pzsandlerextPV.h"
#include "TPZPlasticStepPV.h"
#include "TPZPlasticStep.h"
#include "TPZYCSandlerDimaggio.h"
#include "TPZSandlerDimaggio.h"

void UnaxialLoadingSD();
void ProportionalLoading();
int VerifyTangentSandlerPV();

/*
TPZFNMatrix <6> FromMatToVoight(TPZFNMatrix <9> mat)
{
	TPZFNMatrix <6> voi(6,1,0.);
	int k = 0;
	for (int i = 0; i < 3; i++) {
		for (int j = i ; j < 3; j++) {
			voi(k++,0) = mat(i,j);
		}
	}
	return voi;	
}
*/

void compareplasticsteps()
{
    STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    
    
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
    TPZSandlerExtended SDPV(A,B,C,D,K,G,W,R,N,phi,psi);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    PlasticStepPV.fER =ER;
    
    ofstream outfilePV("comparingPV.txt");
    ofstream outfileErick("comparingErickInitialGuess.txt");
    TPZTensor<STATE> epst,sigma1,sigma2,deps;
    TPZFNMatrix<36> Dep(6,6,0.);
    STATE deltaeps = -0.005;
    
    for(int i=0;i<130;i++)
    {

        PlasticStepPV.ApplyStrainComputeSigma(epst, sigma1);
        outfilePV << -epst.XX()<< " " << -sigma1.XX() << "\n";
        PlasticStepErick.ApplyStrainComputeSigma(epst,sigma2);
        outfileErick << -epst.XX()<< " " << -sigma2.XX() << "\n";
        
        if(i==12)
        {
            deltaeps=0.002;
        }
        if(i==18)
        {
            deltaeps=-0.002;
        }
        if(i==45)
        {
            deltaeps=0.002;
        }
        if(i==64)
        {
            deltaeps=-0.001;
        }
        if (i==130)
        {
            deltaeps=0.001;
        }

        epst.XX()+=deltaeps;
    
        
    }
    
}


void UniaxialSandstone()
{
    
    
    TPZSandlerExtended SDPV;
    SDPV.ReservoirSandstone(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    STATE E= SDPV.fE,nu=SDPV.fnu;
    ER.SetUp(E,nu);
    PlasticStepPV.fER =ER;
    
    ofstream outfilePV("sandstone.txt");
    TPZTensor<STATE> epst,sigma,deps;

    STATE deltaeps = -0.15/80;
    
    for(int i=0;i<140;i++)
    {
        
        PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
        outfilePV << -epst.XX()<< " " << -sigma.XX() << "\n";
        
        if(i==70)
        {
            deltaeps=0.075/80;
        }
     /*  if(i==87)
        {
            deltaeps=-0.15/80;
        }
      */
        if(i==118)
        {
            deltaeps=-0.15/80;
        }
        


        epst.XX()+=deltaeps;
        
        
    }
    
}


void VolumetricTest()
{
    
    
    TPZSandlerExtended SDPV;
    //SDPV.ReservoirSandstone(SDPV);
    SDPV.MCormicRanchSand(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    STATE E= SDPV.fE,nu=SDPV.fnu;
    ER.SetUp(E,nu);
    PlasticStepPV.fER =ER;
    
    ofstream outfilePV("sandstonevol.txt");
    TPZTensor<STATE> epst,sigma,deps;
    
    STATE deltaeps = -0.4/100;
    
    for(int i=0;i<140;i++)
    {
        
        PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
        outfilePV << -epst.I1()<< " " << -sigma.I1() << "\n";
        
         if(i==140)
        {
            deltaeps=0.4/100;
        }
        /* if(i==87)
         {
         deltaeps=-0.15/80;
         }
         
        if(i==118)
        {
            deltaeps=-0.15/80;
        }
        
        */
        
        epst.XX()+=deltaeps;
        epst.YY()+=deltaeps;
        epst.ZZ()+=deltaeps;
        
        
    }
    
}

void ErickTaylorCheck(TPZTensor<REAL> eps, TPZTensor<REAL> deps)
{
	STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
	STATE G=E/(2.*(1.+nu));
	STATE K=E/(3.*(1.-2*nu));
	TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
  PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
	TPZTensor<STATE> sigma,eps1,eps2, Sigma1,Sigma2;
	TPZStack <REAL> coef;
	
	TPZFNMatrix<36> Dep(6,6,0.);
	
	TPZPlasticState<REAL> state;
	state = PlasticStepErick.GetState();
	PlasticStepErick.ApplyStrainComputeDep(eps,sigma,Dep);
	sigma.Print(cout);
	//Dep.Print("DepErick:");
	PlasticStepErick.SetState(state);
	
	REAL scale = 0.01;
	REAL alphatable[] = {0.1,0.2,0.3,0.4,0.5,0.6};
	for (int i = 0; i < 6; i++) {
		alphatable[i] *= scale;
	}
	for (int ia = 0 ; ia < 5; ia++) {
		REAL alpha1 = alphatable[ia];
		REAL alpha2 = alphatable[ia+1];
		eps1.Scale(0.);
		eps2.Scale(0.);
		eps1 = eps;
		eps2 = eps;
		eps1.Add(deps, alpha1);
		eps2.Add(deps, alpha2);
				
		PlasticStepErick.ApplyStrainComputeSigma(eps1,Sigma1);
		PlasticStepErick.SetState(state);
		
		PlasticStepErick.ApplyStrainComputeSigma(eps2,Sigma2);
		PlasticStepErick.SetState(state);
		
		TPZFNMatrix <9> depsMat(3,3,0.);
		TPZFNMatrix <6> deps1(6,1,0.);
		depsMat = deps;
		deps1 = FromMatToVoight(depsMat);
		
		TPZFNMatrix <6> tanmult1(6,1,0.), tanmult2(6,1,0.);
		Dep.Multiply(deps1, tanmult1);
		Dep.Multiply(deps1, tanmult2);
		
		for (int i = 0 ; i < 6; i++) {
			tanmult1(i,0) *= alpha1;
			tanmult2(i,0) *= alpha2;
		}
		
		TPZFNMatrix <9> SigMatTemp33(3,3,0.);
		TPZFNMatrix <6> sigprMat(6,1,0.),sigpr1Mat(6,1,0.),sigpr2Mat(6,1,0.);
		SigMatTemp33 = sigma;
		sigprMat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma1;
		sigpr1Mat = FromMatToVoight(SigMatTemp33);
		SigMatTemp33 = Sigma2;
		sigpr2Mat = FromMatToVoight(SigMatTemp33);
		
		TPZFNMatrix<6> error1(6,1,0.), error2(6,1,0.);
		for (int i = 0 ; i < 6; i++) {
			error1(i,0) = sigpr1Mat(i,0) - sigprMat(i,0) - tanmult1(i,0);
			error2(i,0) = sigpr2Mat(i,0) - sigprMat(i,0) - tanmult2(i,0);
		}
		
		REAL n;
		REAL norm1, norm2;
		norm1 = NormVecOfMat(error1);
		norm2 = NormVecOfMat(error2);
		n = ( log(norm1) - log(norm2) ) / ( log(alpha1) - log(alpha2) );
		coef.push_back(n);
		
	}
	std::cout << "coef = " << coef << std::endl;
}

void comparingDep()
{
    STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    TPZFNMatrix<36> Celast(6,6,0.),Dep(6,6,0.);
    Celast(0,0)=K+(4./3.)*G;Celast(0,3)=K-(2./3.)*G;Celast(0,5)=K-(2./3.)*G;
    Celast(3,0)=K-(2./3.)*G;Celast(3,3)=K+(4./3.)*G;Celast(3,5)=K-(2./3.)*G;
    Celast(5,0)=K-(2./3.)*G;Celast(5,3)=K-(2./3.)*G;Celast(5,5)=K+(4./3.)*G;
    Celast(1,1)=2*G;
    Celast(2,2)=2*G;
    Celast(4,4)=2*G;
    cout << "\n Elastic Contitutive Matrix  C  =  "<< Celast <<endl;
    TPZTensor<STATE> eps, deps,sigma1,sigma2;
    
		eps.XX() = -0.0001;
		eps.YY() = -0.0002;
		eps.ZZ() = -0.0003;
//    eps.XY() = -0.001;
//    eps.XZ()=-0.001;
//    eps.YZ()=-0.001;
    

    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
    TPZSandlerExtended SDPV(A,B,C,D,K,G,W,R,N,phi,psi);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    PlasticStepPV.fER = ER;
    
    cout << "\n C elastic" << Celast <<endl;
#ifdef LOG4CXX
{
    std::stringstream outfile;
    outfile<<"\n Comparing Dep Matrix "<<std::endl;
    
    LOGPZ_DEBUG(logger,outfile.str());
    
}
#endif
    
    PlasticStepPV.ApplyStrainComputeDep(eps, sigma1,Dep);
    cout << "\n Dep PV" << Dep <<endl;
		cout << "\nSigmaPV:" << endl;
		sigma1.Print(cout);

    PlasticStepErick.ApplyStrainComputeDep(eps,sigma2,Dep);
    cout << "\n Dep Erick" << Dep <<endl;
		cout << "\nSigmaErick:" << endl;
	
		sigma2.Print(cout);
		sigma1-=sigma2;
		std::cout << "\nErro entre Sigmas:" << endl;
		sigma1.Print(cout);
	
	/// TAYLOR CHECKS
		deps.XX() = -0.00001;
//		deps.YY() = -0.00001;
//		deps.ZZ() = -0.00001;
//		deps.XY() = -0.00001;
//		deps.XZ() = -0.00002;
//		deps.YZ() = -0.00003;
//    void TPZSandlerExtended::TaylorCheckProjectSigma(const TPZVec<STATE> &sigtrial, STATE kprev, TPZVec<STATE> &xnorm, TPZVec<STATE> &errnorm) const
    {
        TPZManVector<STATE,3> sigtrial(3);
        sigtrial[0] = -0.032;
        sigtrial[1] = -0.04;
        sigtrial[2] = -0.048;
        
        STATE kprev = k0;
        TPZManVector<REAL,11> xnorm,errnorm;
        SDPV.TaylorCheckProjectSigma(sigtrial, kprev, xnorm, errnorm);
        std::cout << xnorm << std::endl;
        std::cout << errnorm << std::endl;
	}
		// CheckConv PV
		PlasticStepPV.TaylorCheck(eps, deps, k0);
		
		// CheckConv do Erick
		ErickTaylorCheck(eps,deps);
}

void SurfacePlot()
{
    TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3),yield(2);
    STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    TPZSandlerExtended materialmodel(A, B,C, D,K,G,W,R, N,phi,psi);
    
    ofstream outfile("surfaceplot.txt");
        ofstream outfile2("surfaceplot2.txt");
    
    deleps[0]= -0.005;
    
    for (int k=0; k<3;k++) {
        deltaepsP[k]=0.;
        epsT[0]=0.;
        epsPnext[0]=0.;
        epssol[0]=0.;
    }
    
    STATE kproj,kprev,epv=0.,I1,J2;
    TPZFMatrix<STATE> GradSigma;
    STATE temp1,temp3;
    
    kproj=0.;
    kprev=0.13304;
    //materialmodel.Firstk(epv,kprev);
    for(int i=0;i<19;i++)
    {
        
        for (int k=0; k<3;k++) {
            epssol[k]=epsT[k]-epsPnext[k];
        }
        
        materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial);
        //materialmodel.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
        materialmodel.ProjectSigma(sigtrial,kprev,sigproj,kproj);
        outfile << -epsT[0]<< " " << -sigproj[0] << "\n";
        //materialmodel.YieldFunction(sigproj, kproj, yield);
        materialmodel.ComputeI1(sigtrial,I1);
        materialmodel.ComputeJ2(sigtrial,J2);
        
         temp1=(-I1+kprev)/(-R*(A-C*exp(B*I1)));
         temp3=(sqrt(J2))/((A-C*exp(B*I1)));
        
        outfile << -I1<<  " " << (A-C*exp(B*I1)) << "\n";
        outfile2<< -I1<<  " "  << temp1*temp1+temp3*temp3 << "\n";
        
        if(i==12)
        {
            deleps[0]=-0.002;
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
    

    
    
    
}


void TwoLoadings()
{
    STATE E=100,nu=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;
    STATE G=E/(2.*(1.+nu));
    STATE K=E/(3.*(1.-2*nu));
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
    TPZSandlerExtended SDPV(A,B,C,D,K,G,W,R,N,phi,psi);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    PlasticStepPV.fER =ER;
    
    ofstream outfilePV("comparingPVFossumTest.txt");
    ofstream outfileErick("comparingErickFossumTest.txt");
    
    TPZTensor<STATE> epst,sigma,deps;
    TPZFNMatrix<36> Dep(6,6,0.);
    STATE deltaeps = -0.1/40.;
    
    for(int i=0;i<80;i++)
    {
        

        
       /*
        
        if(i==40)
        {
            deltaeps = 0.04/40.;
        }
        
        */
        if(i<40)
        {
            PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
            outfilePV << -epst.XX()<< " " << -sigma.XX() << "\n";
            //PlasticStepErick.ApplyStrainComputeSigma(epst,sigma);
            //outfileErick << -epst.XX()<< " " << -sigma.XX() << "\n";
            epst.XX()+=deltaeps;
        }
        else
        {
            PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
            outfileErick << -epst.XY()<< " " << -sigma.XY() << "\n";
            outfilePV << -epst.XX()<< " " << -sigma.XX() << "\n";
            //PlasticStepErick.ApplyStrainComputeSigma(epst,sigma);
            //outfileErick << -epst.XY()<< " " << -sigma.XY() << "\n";
            epst.XY()+=deltaeps;
        }
        
        
        
    }
    
}



int main()
{
    InitializePZLOG();
    {
        std::ofstream out("cadeamerda.txt");
    }
    VerifyTangentSandlerPV();
//    UniaxialSandstone();
    //VolumetricTest();
   // TwoLoadings();
    //SurfacePlot();
//    compareplasticsteps();
    comparingDep();
    
 /*   TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.SetUp(0.25, 100, 0.25, 0.67, 0.18, 2.5, 0.67, 0.066);
    //PlasticStepErick.McCormicRanchSand(PlasticStepErick);
//    TPZSandlerExtended SDPV(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
//    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(0.13045);
//    PlasticStepPV.fYC = SDPV;
//    TPZElasticResponse ER;
//    ER.SetUp(100., 0.25);
//    PlasticStepPV.fER =ER;
    
    TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3);
    
    deleps[0]= -0.005;
    
    for (int k=0; k<3;k++) {
        deltaepsP[k]=0.;
        epsT[0]=0.;
        epsPnext[0]=0.;
        epssol[0]=0.;
    }
    
    STATE kproj,kprev,epv=0.;
    TPZFMatrix<STATE> GradSigma;
    
    ofstream outfile("uniaxialshearPV.txt");
    ofstream outfile2("uniaxialshearPV2.txt");
    ofstream outfile3("I1sqrtJ2a.txt");
    ofstream outfile4("I1sqrtJ2b.txt");
    ofstream outfile5("surface1.txt");
    ofstream outfile6("surface2.txt");
    TPZTensor<STATE> epst,sigma,deps;
    TPZFNMatrix<36> Dep(6,6,0.);
    
    STATE deltaeps = -0.04/40;
    epst.XX()=0;
    epst.XY()=0;
    epst.XZ()=0;
    epst.YY()=0;
    epst.YZ()=0;
    epst.ZZ()=0;
    
    kproj=0.;
    kprev=0.13304;
    //materialmodel.Firstk(epv,kprev);
    
    STATE K=66.6666666666666666666666666666667;
    STATE G=40.;
    TPZFNMatrix<36> Celast(6,6,0.);
    Celast(0,0)=K+(4./3.)*G;Celast(0,3)=K-(2./3.)*G;Celast(0,5)=K-(2./3.)*G;
    Celast(3,0)=K-(2./3.)*G;Celast(3,3)=K+(4./3.)*G;Celast(3,5)=K-(2./3.)*G;
    Celast(5,0)=K-(2./3.)*G;Celast(5,3)=K-(2./3.)*G;Celast(5,5)=K+(4./3.)*G;
    Celast(1,1)=G;
    Celast(2,2)=G;
    Celast(4,4)=G;
    cout << "\n Elastic Contitutive Matrix  C  =  "<< Celast <<endl;
    //PlasticStepErick.SetUp(<#REAL poisson#>, <#REAL E#>, <#REAL A#>, <#REAL B#>, <#REAL C#>, <#REAL R#>, <#REAL D#>, <#REAL W#>);
    PlasticStepErick.SetUp(0.2524,22547.0,689.2, 3.94e-4, 675.2, 28.0, 1.47e-3, 0.08);
    REAL A=689.2,B=3.94e-4,D=1.47e-3,C=675.2,E=22547,nu=0.2524,R=28.,W=0.08;
    G=E/(2.*(1.+nu));
    K=E/(3.*(1.-2*nu));
//TPZSandlerExtended(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi):
    TPZSandlerExtended SDPV(A,B,C, D,K,G,W,R, 0,6,1);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(-8.05);
    PlasticStepPV.fYC = SDPV;
    TPZElasticResponse ER;
    ER.SetUp(E, nu);
    PlasticStepPV.fER =ER;
    STATE epspp=0.,kk;
    
    TPZVec<STATE> strain(3),loading(3),tyild(2);
    strain[0]=-0.04/40;
    strain[1]=0.;
    strain[2]=0.;
    
    SDPV.Firstk(epspp, kk);
    
    
    for(int i=0;i<80;i++)
    {
        
        //PlasticStepPV.ApplyStrainComputeDep(epst, sigma, Dep);
        PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
        cout<<"\n i "<< i;
        //PlasticStepErick.ApplyStrainComputeDep(epst, sigma, Dep);
        //PlasticStepErick.ApplyStrainComputeSigma(epst, sigma);
        

        
        if(i>=40)
        {
            outfile2 << epst.XY()<< " " << sigma.XY() << "\n";
            outfile << -epst.XX()<< " " << -sigma.XX() << "\n";
            outfile3 << -sigma.I1() << " "<<sqrt(sigma.J2())<<"\n";
            outfile4 << -sigma.I1() << " "<<sqrt(sigma.J2())-(A-C*exp(B*-sigma.I1()))<<"\n";
        }
        else
        {
            outfile << -epst.XX()<< " " << -sigma.XX() << "\n";
            outfile3 << -sigma.I1() << " "<<sqrt(sigma.J2())<<"\n";
            outfile4 << -sigma.I1() << " "<<sqrt(sigma.J2())-(A-C*exp(B*-sigma.I1()))<<"\n";
        }
    
        if(i==40)
        {
            deltaeps=-0.04/40;
            deltaeps*=-1;
        }
        
        if(i>=40)
        {
            epst.XY()+=deltaeps;
        }
        else
        {
           epst.XX()+=deltaeps;
        }
        
    }
    
 
*/

//    for(int i=0;i<15;i++)
//    {
//        
//        //PlasticStepPV.ApplyStrainComputeDep(epsT, sigma, Dep);
//        PlasticStepPV.ApplyStrainComputeSigma(epsT, sigma);
//        
//        outfile << -epsT[0]<< " " << -sigma[0] << "\n";
//        //cout << "\n Dep Diogo = "<< Dep <<endl;
//        
//        //PlasticStepErick.ApplyStrainComputeDep(epsT, sigma, Dep);
//        PlasticStepErick.ApplyStrainComputeSigma(epsT, sigma);
//        
//        outfile2 << -epsT[0]<< " " << -sigma[0] << "\n";
//        //cout << "\n Dep Erick  = "<< Dep <<endl;
//        
//        if(i==12)
//        {
//            deltaeps*=-1;
//        }
//        epsT.XX()+=deltaeps;
//    }

    
    
    
    
  //  UnaxialLoadingSD();
    
    return 0;
}

int VerifyTangentSandlerPV()
    {
    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
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
    
    TPZManVector<STATE,3> deltaeps(3),eps(3),sigma(3),deltasigma(3);
    

    deltaeps[0]= -0.00135;
    deltaeps[1]=0;
    deltaeps[2]=0;
    eps=deltaeps;
    STATE kprev=0.,epspv1=0.;
    materialmodel.Firstk(epspv1,kprev);
    TPZManVector<STATE,3> Tensor2(3);
    Tensor2[0]=-0.032;
    Tensor2[1]=-0.04;
    Tensor2[2]=-0.048;
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
    
    cout << "Teste de D(DistF2)/D(theta,beta,k)\n";
        {
            STATE kproj;
            TPZManVector<STATE,3> sigproj(3);
            materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
            STATE theta,beta;
            materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);
            materialmodel.TaylorCheckDistF2(Tensor2, theta, beta, kproj, kprev, xnorm, errnorm);
        }
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
        cout << "Teste da derivada D(ResF2)/D(theta,beta,k)\n";
    materialmodel.TaylorCheckDDistF2(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
    cout << "Teste da derivada D(ResF2)/D(sigtrial)\n";
    {
        STATE kproj;
        TPZManVector<STATE,3> sigproj(3);
        materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
        STATE theta,beta;
        materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);

        materialmodel.TaylorCheckDDistF2DSigtrial(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
    }
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "xnorm " << xnorm << endl;
    cout << "errnorm " << errnorm << endl;
    cout << "convergence rate " << converge << endl;
    
    materialmodel.TaylorCheckDF1Cart(xi, beta, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
        
    
        cout << "Teste de D(F2cart)/D(theta,beta,k)\n";
        {
            STATE kproj;
            TPZManVector<STATE,3> sigproj(3);
            materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
            STATE theta,beta;
            materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);
            materialmodel.TaylorCheckDF2Cart(theta, beta, kproj, xnorm, errnorm);
        }
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "\n xnorm " << xnorm << endl;
    cout << "\n errnorm " << errnorm << endl;
    cout << "\n convergence rate " << converge << endl;
    
        std::cout << "Teste de D(sigproj)/D(sigtrial) geral\n";
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
    
        std::cout << "Teste de D(sigproj)/D(sigtrial) para funcao F2\n";
        materialmodel.TaylorCheckProjectF2(Tensor2, kprev, xnorm, errnorm);
        materialmodel.ConvergenceRate(xnorm, errnorm, converge);
        cout << "\n xnorm " << xnorm << endl;
        cout << "\n errnorm " << errnorm << endl;
        cout << "\n convergence rate " << converge << endl;
        
        std::cout << "Teste de D(theta,beta,k)/D(sigtrial) para funcao F2\n";
    materialmodel.TaylorCheckDtbkDsigtrial(Tensor2, kprev, xnorm, errnorm);
    materialmodel.ConvergenceRate(xnorm, errnorm, converge);
    cout << "xnorm " << xnorm << endl;
    cout << "errnorm " << errnorm << endl;
    cout << "convergence rate " << converge << endl;
        
   
	return 0;
}

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



//UnaxialLoadingSD();
//ProportionalLoading();



/*
 TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3);
 
 //   TPZSandlerExtended::TPZSandlerExtended(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi);
 REAL E = 22547.,nu = 0.2524;
 STATE K = E/(3.*(1.-2.*nu));
 STATE G = E/(2.*(1.+nu));
 TPZSandlerExtended materialmodel2(689.2, 3.94e-4,675.2,1.47e-3,K,G,0.08,28., 0,6.,1.);
 
 ofstream outfile2("fossum1.txt");
 ofstream outfile3("fossumI1sqrtJ2.txt");
 
 deleps[0]= -0.04/40;
 
 for (int k=0; k<3;k++) {
 deltaepsP[k]=0.;
 epsT[0]=0.;
 epsPnext[0]=0.;
 epssol[0]=0.;
 }
 
 STATE kproj,kprev,epv=0.,I1,sqrtJ2;
 TPZFMatrix<STATE> GradSigma;
 
 kproj=0.;
 kprev=-8.05;
 //materialmozdel.Firstk(epv,kprev);
 for(int i=0;i<=40;i++)
 {
 
 for (int k=0; k<3;k++) {
 epssol[k]=epsT[k]-epsPnext[k];
 }
 
 materialmodel2.ApplyStrainComputeElasticStress(epssol, sigtrial);
 materialmodel2.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
 materialmodel2.ProjectSigma(sigtrial,kprev,sigproj,kproj);
 outfile2 << -epsT[0]<< " " << -sigproj[0] << "\n";
 materialmodel2.ComputeI1(sigproj,I1);
 materialmodel2.ComputeJ2(sigproj,sqrtJ2);
 sqrtJ2=sqrt(sqrtJ2);
 outfile3 << I1<< " " << sqrtJ2 << "\n";
 
 if(i==40)
 {
 deleps[0]=0.04/40;
 deleps[1]=0.04/40;
 }
 
 for (int k=0; k<3;k++) {
 deltasigma[k]=sigtrial[k]-sigproj[k];
 }
 
 materialmodel2.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);
 
 for (int k=0; k<3;k++) {
 epsPnext[k]+=deltaepsP[k];
 epsT[k]+=deleps[k];
 }
 kprev=kproj;
 
 }
 */

void UnaxialLoadingSD()
{
    TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3);
    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
  
    ofstream outfile("fortransduniaxial.txt");
    
    deleps[0]= -0.005;
    
    for (int k=0; k<3;k++) {
        deltaepsP[k]=0.;
        epsT[0]=0.;
        epsPnext[0]=0.;
        epssol[0]=0.;
    }
    
    STATE kproj,kprev,epv=0.;
    TPZFMatrix<STATE> GradSigma;
    
    kproj=0.;
    kprev=0.13304;
    //materialmodel.Firstk(epv,kprev);
    for(int i=0;i<19;i++)
    {
        
        for (int k=0; k<3;k++) {
            epssol[k]=epsT[k]-epsPnext[k];
        }
        
        materialmodel.ApplyStrainComputeElasticStress(epssol, sigtrial);
        materialmodel.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
        materialmodel.ProjectSigma(sigtrial,kprev,sigproj,kproj);
        outfile << -epsT[0]<< " " << -sigproj[0] << "\n";
        
        if(i==12)
        {
            deleps[0]=-0.002;
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
    
}

void ProportionalLoading()
{
 
    TPZManVector<STATE,3> epsPnext(3),epsT(3),deleps(3),epssol(3),deltaepsP(3),sigproj(3),sigtrial(3),deltasigma(3),dsig(3);
    TPZSandlerExtended materialmodel(0.25, 0.67,0.18, 0.67,66.67,40.,0.066,2.5, 0,0,1);
    
    ofstream outfile("fortransdproportional04.txt");
    
    dsig[0]= -1.4/11.;
    dsig[1]= dsig[0];
    dsig[2]= dsig[0];
    
    materialmodel.ApplyStressComputeElasticStrain(dsig,deleps);
    
    for (int k=0; k<3;k++) {
        deltaepsP[k]=0.;
        epsT[0]=0.;
        epsPnext[0]=0.;
        epssol[0]=0.;
    }
    
    STATE kproj,kprev,epv=0.;
    TPZFMatrix<STATE> GradSigma;
    
    kproj=0.;
    //kprev=0.13304;
    REAL epsp=0;
    materialmodel.Firstk(epsp,kprev);
    for(int i=0;i<11;i++)
    {
        
        for (int k=0; k<3;k++) {
            epssol[k]=epsT[k]-epsPnext[k];
        }
        
        materialmodel.ApplyStrainComputeElasticStress(epsT, sigtrial);
        materialmodel.ProjectSigmaDep(sigtrial,kprev,sigproj,kproj,GradSigma);
        materialmodel.ProjectSigma(sigtrial,kprev,sigproj,kproj);
        outfile << -epsT[0]<< " " << -sigproj[0] << "\n";
        
        
        for (int k=0; k<3;k++) {
            deltasigma[k]=sigtrial[k]-sigproj[k];
        }
        
        materialmodel.ApplyStressComputeElasticStrain(deltasigma, deltaepsP);
        
        for (int k=0; k<3;k++) {
            epsPnext[k]+=deltaepsP[k];
            
        }
        dsig[0]+= -1.4/11.;
        //dsig[1]= 0.4*dsig[0];
        //dsig[2]= 0.4*dsig[0];
        materialmodel.ApplyStressComputeElasticStrain(dsig,deleps);
        
        epsT[0]=deleps[0];
        epsT[1]=deleps[1];
        epsT[2]=deleps[2];
        kprev=kproj;
        
    }
    

}

