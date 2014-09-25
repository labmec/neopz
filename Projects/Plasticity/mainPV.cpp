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
#include "TPZTimer.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

#ifdef LOG4CXX
static LoggerPtr loggerConvTest(Logger::getLogger("ConvTest"));
#endif

using namespace pzshape; // needed for TPZShapeCube and related classes

#include <math.h>

//#define MACOS
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

#include "TPZYCMohrCoulombPV.h"

void UnaxialLoadingSD();
void ProportionalLoading();
int VerifyTangentSandlerPV();
void ErickTaylorCheck(TPZTensor<REAL> eps, TPZTensor<REAL> deps);
void CheckDepConv();
void UniaxialLoadingPV(TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> PlasticStepPV);
void UniaxialLoadingPV(TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV,TPZVec<TPZPlasticState<REAL> > &state);
void UniaxialLoadingErick(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick,TPZVec<TPZPlasticState<REAL> > &state);
void I1vsSqrtJ2();
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
    PlasticStepPV.fER =SDPV.GetElasticResponse();
    
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


void I1vsSqrtJ2()
{
  
    
    TPZSandlerExtended SDPV;
    SDPV.MCormicRanchSand(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    PlasticStepPV.fER =SDPV.GetElasticResponse();
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.McCormicRanchSand(PlasticStepErick);
    
    ofstream outfilePV("i1vsj22.txt");
    ofstream outfilePVtrial("uniaxialsolu.txt");
    ofstream outfileFf("Ff(I1,0).txt");
    ofstream myfile("teste.txt");

    TPZTensor<STATE> epst,sigma,sigma2,deps,stresstrial;

    STATE deltaeps = -0.005;
    STATE sqrtj2,f2;

    for(int i=0;i<20;i++)
    {        
        
        //ER.Compute(epst, stresstrial);

        //PlasticStepErick.ApplyStrainComputeSigma(epst,sigma2);
        //outfilePVtrial << sigma2.I1() << " " << sqrt(sigma2.J2()) <<"\n";
        
        
        PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
        outfilePV << -sigma.I1()<< " " << sqrt(sigma.J2()) << "\n";
        outfilePVtrial << -epst.XX() << " " << -sigma.XX() <<"\n";
      
        
        k0 = PlasticStepPV.fN.Alpha();
        STATE i1=k0,X,tmp;
        sqrtj2=sqrt(sigma.J2());
        tmp=1.;
        TPZTensor<STATE>::TPZDecomposed principal;
        sigma.EigenSystem(principal);
        TPZVec<STATE> pt(3,0.),yield(2,0.);
        pt=principal.fEigenvalues;
        SDPV.YieldFunction(pt, k0, yield);
        if(yield[1]<1.e-6)
        {
            for(int j =0;j<100;j++)
            {
                X=SDPV.GetX(k0);
                STATE val=(X-k0)*(X-k0)-(i1-k0)*(i1-k0);
                if (val<0) {
                val=0;
                }
                f2=(1./SDPV.GetR())*sqrt(val);

                if(tmp!=0)
                {
                    myfile << -i1 << " " <<f2 << "\n";
                }
    
                tmp=f2;
                i1+=SDPV.GetX(k0)/60.;
            }
            myfile << "\n\n\n";
        }
        
        i1=-3;
        for(int ii=0;ii<30;ii++)
        {
            outfileFf << -i1 << " "<< SDPV.GetF(i1) <<"\n";
            i1+=0.1205;
            
        }
        outfileFf << "\n\n";
        
         if(i==12)
        {
            deltaeps=0.002;
        }
//        if(i==15)
//        {
//            deltaeps=0.001;
//        }
//        if(i==16)
//        {
//            deltaeps=0.002;
//        }
        
        epst.XX()+=deltaeps;

        
    }
 
}


void UniaxialLoadingApplyStrainComputeDep()
{
    
    
    TPZSandlerExtended SDPV;
    SDPV.MCormicRanchSand(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;

    PlasticStepPV.fER =SDPV.GetElasticResponse();
    
    
    STATE E =100,nu=0.25;
    
    
    STATE A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,N=0.,phi=0,psi=1.0;

    
   
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
    //TPZSandlerExtended SDPV(A,B,C,D,K,G,W,R,N,phi,psi);
    
    //TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    //PlasticStepErick.McCormicRanchSand(PlasticStepErick);
    
    ofstream uniaxialPV("XuniaxialPV.txt");
    ofstream uniaxialER("XuniaxialERICK.txt");
    
    TPZTensor<STATE> epst,sigma1,sigma2,deps,stresstrial;
    
    STATE deltaeps = -0.005;
    TPZFNMatrix<36> DepPV(6,6,0.),DepER(6,6,0.);
    STATE DsigxDepsx,tempsig1=0,tempsig2=0,tempeps=0;
    for(int i=0;i<20;i++)
    {

       PlasticStepErick.ApplyStrainComputeDep(epst,sigma1,DepER);
       uniaxialER << -epst.XX() << " " << -sigma1.XX() <<"\n";
        
       PlasticStepPV.ApplyStrainComputeDep(epst, sigma2,DepPV);
       uniaxialPV << -epst.XX() << " " << -sigma2.XX() << "\n";
        
#ifdef LOG4CXX
        {
            std::stringstream outfile;
            
            DsigxDepsx = (sigma1.XX()-tempsig1)/(epst.XX()-tempeps);
            outfile<<"\n DsigxDepsx Erick "<< DsigxDepsx <<std::endl;
            outfile<<"\n Dep ERICK \n "<< DepER <<std::endl;
            
            DsigxDepsx = (sigma2.XX()-tempsig2)/(epst.XX()-tempeps);
            outfile<<"\n DsigxDepsx PV "<< DsigxDepsx <<std::endl;
            outfile<<"\n Dep PV \n "<< DepPV <<std::endl;
            LOGPZ_DEBUG(logger,outfile.str());
            
        }
#endif
        tempsig1=sigma1.XX();
        tempsig2=sigma2.XX();
        tempeps=epst.XX();
        if(i==12)
        {
            deltaeps=0.002;
        }
        deps.XX()=deltaeps;
        epst.XX()+=deltaeps;
				TPZManVector<REAL,3> conv;
        PlasticStepPV.TaylorCheck(epst, deps, k0, conv);
		
		// CheckConv do Erick
		ErickTaylorCheck(epst,deps);
        
    }
    
}




void ReadData(const std::string &filename)
{
    std::ifstream input(filename.c_str());
    if (!input) {
        DebugStop();
    }
    ofstream LabUCS("LABDATA_EpsXvsSigX.txt");
    //ofstream LabRR("LABDATA_EpsRvsSigR.txt");
    ofstream LabVol("LABDATA_EpsVvsSigV.txt");
    ofstream TempoVsEpsV("LABDATA_TempoVsEpsV.txt");
    ofstream TempoVsEpsAxi("LABDATA_TempoVsEpsAxi.txt");
    ofstream TempoVsEpsLat("LABDATA_TempoVsEpsLat.txt");
    ofstream TempoVsSigV("LABDATA_TempoVsSigV.txt");
    ofstream TempoVsSigAxi("LABDATA_TempoVsSigAxi.txt");
    ofstream TempoVsSigLat("LABDATA_TempoVsSigLat.txt");
    ofstream I1vsSqtJ2("LABDATA_I1vsSqtJ2.txt");
     ofstream EpsLatVsSigLat("LABDATA_EpsLatVsSigLat.txt");
    ofstream I1vsEPV("LABDATA_I1vsEPV.txt");
    int numlines = 0;
    char buf[1024];
    input.getline(buf , 1024);
    //vector<double> X,sigx,time,sigdev,sigr,epsx,epsr,epsv;
    STATE X , Sig_Axial,Sig_Axial_Desv,Sig_Conf,Def_Axial,Def_Lateral,Tempo,Def_Volume;
    while (input) {
       // input >> Tempo1 >> Sig_Axial_Desv >> Tempo1 >> Def_Axial >> Tempo1 >> Def_Lateral >> Tempo1 >> Def_Volume >> Tempo1 >> Def_Volume;
//        X	       Sig Axial Total	  Tempo	       Sig Axial Desv	Tempo	        Sig Conf	Tempo	       Def Axial	Tempo	      Def Lateral	Tempo
        input >> X>>Sig_Axial>>Tempo>>Sig_Axial_Desv>>Tempo>>Sig_Conf>>Tempo>>Def_Axial>>Tempo>>Def_Lateral>>Tempo>>Def_Volume;
        if (!input) {
            break;
        }
        
        //LabUCS << Def_Axial/10. << " " << Sig_Axial_Desv <<"\n";
        //LabRR << eps_r << " " << sig_r << "\n";
        LabVol << /*(Def_Axial+2*Def_Lateral)*/(Def_Axial+2*Def_Lateral) << " " << (Sig_Axial+2*Sig_Conf) << "\n";
        LabUCS << /*(Def_Axial+2*Def_Lateral)*/Def_Axial << " " << Sig_Axial << "\n";
        TempoVsEpsV << Tempo<< " " << (Def_Axial+2*Def_Lateral) << "\n";
        TempoVsEpsAxi << Tempo << " " << Def_Axial  << "\n";
        TempoVsEpsLat << Tempo << " " << Def_Lateral << "\n";
        TempoVsSigV << Tempo << " " << (Sig_Axial+2*Sig_Conf)  << "\n";
        TempoVsSigAxi << Tempo << " " << Sig_Axial << "\n";
        TempoVsSigLat << Tempo << " " << Sig_Conf << "\n";
        EpsLatVsSigLat << Def_Lateral << " " << Sig_Conf << "\n";
        
        TPZTensor<STATE> stress;
        stress.XX()=Sig_Axial;
        stress.YY()=Sig_Conf;
        stress.ZZ()=Sig_Conf;
        stress.XY()=Sig_Axial_Desv;
        I1vsSqtJ2 << stress.I1() << " " << sqrt(stress.J2()) << "\n";
        cout << stress.I1() << " " << sqrt(stress.J2()) << "\n";
         I1vsEPV << stress.I1() <<" "<<(Def_Axial+2*Def_Lateral)/100.<<"\n";
//       /* X.push_back(x);
//        sigx.push_back(sig_ax_t);
//        time.push_back(tempo);
//        sigdev.push_back(sig_ax_dev);
//        sigr.push_back(sig_r);
//        epsx.push_back(eps_ax);
//        epsr.push_back(eps_r);
//        epsv.push_back(eps_v);
//        
//     */   
    }

    
    TPZSandlerExtended SDPV;
    SDPV.PreSMat(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
    PlasticStepPV.fYC = SDPV;
    PlasticStepPV.fER =SDPV.GetElasticResponse();
    
    ofstream outfilePV("I1vsSqrtJ2Calc.txt");
    ofstream outfileFf("Ff(I1,0)Ps.txt");
    ofstream myfile("testePs.txt");
    
    TPZTensor<STATE> epst,sigma,deps;
    
//    STATE deltaeps = -1.;
    
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
    PlasticStepErick.PRSMatMPa(PlasticStepErick);
   // PlasticStepErick.Print(cout);
    //PlasticStepErick.UncDeepSandResMPa(PlasticStepErick);
  //  STATE sqrtj2,f2;
//    sigma.XX()=-40.;
//    sigma.XY()=-0.52;
//    sigma.XZ()=-0.52;
//    sigma.YZ()=-0.52;
//    sigma.YY()=-40.;
//    sigma.ZZ()=-40.;
//    PlasticStepErick.ApplyLoad(sigma, epst);
//     PlasticStepErick.Print(cout);
//    for(int i=0;i<160;i++)
//    {
//        
//        
//        PlasticStepErick.ApplyLoad(sigma, epst);
//        //PlasticStepPV.ApplyStrainComputeSigma(epst, sigma);
//        outfilePV << -sigma.I1()<< " " << sqrt(sigma.J2()) << "\n";
//        myfile << -epst.XX()<< " " <<-sigma.XX() << "\n";
//        
//        if(i>40)
//        {
//            sigma.XX()+=-1.;
//        }
//        else
//        {
//            sigma.XX()+=-1.;
//            sigma.YY()+=-1.;
//            sigma.ZZ()+=-1.;
//        }
//        
//        //sigma.XX()+=-1.;
//    }
//    
    
    
}


//REAL TPBrStrainStressDataBase::I1(int index)
//{
//    REAL I1 = this->fSig_Ax[index]+2.*this->fSig_Lat[index];
//    return I1;
//}
//
///// retorno o valor de Sq(J2) para o index
//REAL TPBrStrainStressDataBase::SqJ2(int index)
//{
//    REAL I1 = this->I1(index);
//    REAL sigdesv[3] ={fSig_Ax[index]-I1/3.,fSig_Lat[index]-I1/3.,fSig_Lat[index]-I1/3.};
//    REAL J2 = (sigdesv[0]*sigdesv[0]+sigdesv[1]*sigdesv[1]+sigdesv[2]*sigdesv[2])/2.;
//    return sqrt(J2);
//}
//
///// retorna o valor de Eps Volumetrico para o index
//REAL TPBrStrainStressDataBase::Epsv(int index)
//{
//    REAL epsvol = fEps_Ax[index] + 2 * fEps_Lat[index];
//    return epsvol;
//}
//
///// retorna o valor de Sig Volumetrico para o index
//REAL TPBrStrainStressDataBase::Sigv(int index)
//{
//    REAL sigvol = fSig_Ax[index] + 2 * fSig_Ax[index];
//    return sigvol;
//}
//
//// retorna o valor da envoltoria para o index
//REAL TPBrStrainStressDataBase::F1(int index)
//{
//    REAL I1 = this->I1(index);
//    REAL F1 = fA - fC*exp(fB*I1);
//    return F1;
//}



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
		
		for (int i = 0 ; i < 6; i++)
        {
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
#ifdef LOG4CXX
    {
        std::stringstream outfile;
        outfile<<"\n Check Erick Coef \n "<< coef <<std::endl;
        LOGPZ_DEBUG(logger,outfile.str());
        
    }
#endif
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
        //eps.XY() = -0.001;
        //eps.XZ()=-0.001;
        //eps.YZ()=-0.001;
    

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
    //std::stringstream outfile;
    //outfile<<"\n Comparing Dep Matrix "<<std::endl;
    
    //LOGPZ_DEBUG(logger,outfile.str());
    
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
	TPZManVector<REAL,3> conv;
		PlasticStepPV.TaylorCheck(eps, deps, k0,conv);
		
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

    TPZSandlerExtended SDPV;
    SDPV.SalemLimestone(SDPV);
    STATE epsp=0.,k0;
    SDPV.Firstk(epsp, k0);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(-8.05);
    PlasticStepPV.fYC = SDPV;
    PlasticStepPV.fER =SDPV.GetElasticResponse();
    
    ofstream outfilePV("comparingPVFossumTest.txt");
    ofstream outfileSHEAR("comparingPVFossumTesSHEAR.txt");
    
    TPZTensor<STATE> epst,sigma,deps;
    TPZFNMatrix<36> Dep(6,6,0.);
    STATE deltaeps = -0.04/40.;
    
    for(int i=0;i<80;i++)
    {
        

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
            outfilePV << -epst.XX()<< " " << -sigma.XX() << "\n";
            //PlasticStepErick.ApplyStrainComputeSigma(epst,sigma);
            outfileSHEAR << -epst.XY()<< " " << -sigma.XY() << "\n";
            epst.XY()+=deltaeps;
        }
        
        
        
    }
    
}

void DistFunc2TangentTest()
{
//    TPZSandlerExtended SDPV;
//    SDPV.MCormicRanchSand(SDPV);
//    TPZVec<STATE> pt(3);
//    STATE theta=M_PI/3.,beta=M_PI,k0,epsp0=0.,kprev=0.;
//    TPZFMatrix<STATE> d2distf2(3,3,0.),ddistf2(3,1,0.);
//    SDPV.Firstk(epsp0,k0);
//    pt[0]=0.0160555,pt[1]= -0.0560555,pt[2]= -0.06;
//    
//    SDPV.D2DistFunc2new(pt, theta,beta, k0, d2distf2);
//    cout <<"\n D2DistFunc2New "<<d2distf2<<endl;
//    
//    SDPV.D2DistFunc2(pt, theta,beta, k0, d2distf2);
//    cout <<"\n D2DistFunc2Old "<<d2distf2<<endl;
//
//    TPZManVector<STATE> ddistf2vec(3);
//    SDPV.DDistFunc2new(pt,theta,beta,k0,kprev,ddistf2vec);
//    cout <<"\n ddistf2New "<<ddistf2vec<<endl;
//    
//    
//    SDPV.DDistFunc2(pt, theta, beta, k0, kprev, ddistf2vec);
//    cout <<"\n ddistf2Old "<<ddistf2vec<<endl;
//    
}

void InitialLoad(TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV,TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick)
{
    TPZTensor<REAL> initstress,finalstress,fConfinement,eps;
    finalstress.XX()=-98.8298375333333;
    finalstress.YY()=-98.8298375333333;
    finalstress.ZZ()=-98.8298375333333;
    fConfinement.XX()=-45.9;
    fConfinement.YY()=-62.1;
    fConfinement.ZZ()=-48.2;
    
//    finalstress*=1/10.;
//    PlasticStepErick.ApplyLoad(finalstress, eps);
//    PlasticStepPV.ApplyLoad(finalstress, eps);
    
    PrepareInitialMat(PlasticStepPV, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = fConfinement;
    PrepareInitialMat(PlasticStepPV, initstress, finalstress, 10);
    
    cout <<"\n StatePV = " <<PlasticStepPV.GetState()<<endl;
    
    initstress.XX()=0.;
    initstress.YY()=0.;
    initstress.ZZ()=0.;
    finalstress.XX()=-98.8298375333333;
    finalstress.YY()=-98.8298375333333;
    finalstress.ZZ()=-98.8298375333333;
    fConfinement.XX()=-45.9;
    fConfinement.YY()=-62.1;
    fConfinement.ZZ()=-48.2;
    PrepareInitialMat(PlasticStepErick, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = fConfinement;
    PrepareInitialMat(PlasticStepErick, initstress, finalstress, 10);
    
    cout <<"\n StateFull = " <<PlasticStepErick.GetState()<<endl;

}


int main()
{
 /*   //STATE val = atan2(3,1.);
    REAL poisson = 0.203;
    REAL elast = 29269.;
    REAL A = 152.54;
    REAL B = 0.0015489;
    REAL C = 146.29;
    REAL R = 0.91969;
    REAL D = 0.018768;
    REAL W = 0.006605;
    STATE G=elast/(2.*(1.+poisson));
    STATE K=elast/(3.*(1.-2*poisson));
    TPZSandlerExtended SDPV;
    SDPV.SetUp(A,B, C, D, K, G, W, R, 0, 0, 1.);


    TPZManVector<STATE> sigtr(3),sigproj(3);
    STATE kprev=-226.3,kproj;
    sigtr[0]=100.934;
    sigtr[1]=-61.279;
    sigtr[2]=-251.83;
    SDPV.ProjectF1(sigtr, kprev, sigproj,kproj);
    
    
    return 0;
  */
 
    InitializePZLOG();
    
    UnaxialLoadingSD();
    
    return 0;
    
    //CheckDepConv();
    //ReadData("ensaio_UCS_all_columns.txt");
    //ReadData("ensaio_all_columns.txt");
    //DistFunc2TangentTest();
    //VerifyTangentSandlerPV();
    //comparingDep();
    //I1vsSqrtJ2();
    //TwoLoadings();

//    const int nVarsResidual = 7+2;
//    const int nVarsTensor = 6;
//    typedef TFad<nVarsTensor, REAL> TFAD;
//    typedef TFad<nVarsResidual, TFAD> TFAD_FAD;
    // Aprendendo a usar FAD
    
//    TFad<1,REAL> gx, fx;
//    
//    fx.fastAccessDx(0)=1.;
//    std::cout << fx << std::endl;
//    fx.val()=1.;
//    std::cout << fx << std::endl;
//    gx=cos(fx*fx);
//    
//    std::cout << gx << std::endl;
//    TFad<3,REAL>a,b,c;
//    a.fastAccessDx(0)=1.;
//    a.val()=1.;
//    b.fastAccessDx(1)=1.;
//    b.val()=1.;
//    c.fastAccessDx(2)=1.;
//    c.val()=1.;
//    TFad< 3,REAL > f1,f2,f3;
//    f1 = 2.*a*c*c;
//    f2 = 3.*b*b*a;
//    f3 = a*a*b*c;
//    
//    std::cout << f1 << std::endl;
//    std::cout << f2 << std::endl;
//    std::cout << f3 << std::endl;
    
    
    //comparingDep();
    //UniaxialLoadingApplyStrainComputeDep();
    
    TPZTimer time;
    
    TPZSandlerExtended SDPV;
    STATE epsp,k;
    
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV;
    
    TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> PlasticStepPVMohr;

    
    
    REAL poisson = 0.203;
    REAL Elast = 29269.;
    REAL A = 152.54;
    REAL B = 0.0015489;
    REAL C = 146.29;
    REAL R = 0.91969;
    REAL D = 0.018768;
    REAL W = 0.006605;
    //STATE G=Elast/(2.*(1.+poisson));
    //STATE K=Elast/(3.*(1.-2*poisson));
    
    //STATE Elast=100,poisson=0.25,A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066;
	TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
//	PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
    
    PlasticStepErick.SetUp(poisson, Elast , A, B, C, R, D, W);
    STATE G=Elast/(2.*(1.+poisson));
    STATE K=Elast/(3.*(1.-2*poisson));
    STATE phi=0,psi=1.,N=0.;
//  void TPZSandlerExtended::SetUp(STATE A, STATE B,STATE C, STATE D,STATE K,STATE G,STATE W,STATE R,STATE Phi,STATE N,STATE Psi)
    PlasticStepPV.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
    PlasticStepPV.fER.SetUp(Elast,poisson);
    PlasticStepPV.fYC.Firstk(epsp,k);
//    TPZPlasticState<REAL> state;
//    state.fAlpha=-41.0127;
//    PlasticStepPV.SetState(state);
//    
    InitialLoad(PlasticStepPV,PlasticStepErick);
    
    time.reset();
    time.start();
    TPZVec<TPZPlasticState<REAL> > statePV,stateFull;
    for(int i=1;i<2;i++)
    {
        UniaxialLoadingPV(PlasticStepPV,statePV);
    }
    time.stop();
    
    cout << "\n tempo PV " <<time.seconds() << endl;
    //cout <<"\n STATE PV" <<PlasticStepPV.GetState()<<endl;
    time.reset();

	
    time.start();
    for(int i=1;i<2;i++)
    {
        UniaxialLoadingErick(PlasticStepErick,stateFull);
    }
    time.stop();
    cout << "\n tempo erick " <<time.seconds() << endl;
    TPZTensor<REAL> diff;
    REAL alfa;
    for(int i=0;i<statePV.size();i++)
    {
        diff=statePV[i].fEpsT;
        diff-=stateFull[i].fEpsT;
        cout<<"\n epsT diff fEpsT = "<<diff<<endl;
        diff=statePV[i].fEpsP;
        diff-=stateFull[i].fEpsP;
        cout<<"\n epsT diff fEpsP = "<<diff<<endl;
        alfa=statePV[i].fAlpha;
        alfa-=stateFull[i].fAlpha;
        cout<<"\n epsT diff fAlpha = "<<alfa<<endl;
    }
    

    //cout <<"\n STATE FAD" <<PlasticStepPV.GetState()<<endl;
	  time.reset(); 
	    
	
	//UniaxialLoadingPV(PlasticStepPV,PlasticStepErick);
	//UniaxialLoadingErick(PlasticStepErick);
	return 0;
  
}

int VerifyTangentSandlerPV()
{
   
    
	const REAL A = 0.25, B = 0.67, C = 0.18, D = 0.67, K = 66.67, G = 40., W = 0.066, R = 2.5, Phi = 0., N = 0., Psi = 1.;    
	TPZSandlerExtended materialmodel(A, B, C, D, K, G, W, R, Phi, N, Psi);
	TPZPlasticState<REAL> plasticstate;
	
	TPZElasticResponse ER= materialmodel.GetElasticResponse();
	TPZTensor <STATE> eps1, sigma1, deps1, dsig1;
	sigma1.XX()=0.0128547;
	sigma1.YY()=-0.0564273;
	sigma1.ZZ()=-0.0564273;
	ER.ComputeDeformation(sigma1, eps1);
	deps1.YY() += eps1.YY()/200.;
	ER.Compute(deps1, dsig1);
	dsig1.Print(cout);
	
	STATE kprev=0.,epspv1=0.,knext=0.;
	materialmodel.Firstk(epspv1,kprev);
	
	TPZManVector<STATE,3> Tensor2(3,0.),SigProj(3,0.);
	Tensor2[0] = sigma1.XX();
	Tensor2[1] = sigma1.YY();
	Tensor2[2] = sigma1.ZZ();
	//    Tensor2[0]=-0.032;
	//    Tensor2[1]=-0.04;
	//    Tensor2[2]=-0.048;
	materialmodel.CheckCoordinateTransformation(Tensor2);
	
	//Yield antes do project
	TPZManVector<STATE,2> yield(2);
	materialmodel.YieldFunction(Tensor2,kprev,yield);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\n sigtrial = " << Tensor2 << endl;
		str << "\n k = " << kprev << endl;
		str << "\n yield = " << yield << endl; 
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	//Yield depois do project
	materialmodel.ProjectSigma(Tensor2,kprev,SigProj,knext);
	materialmodel.YieldFunction(SigProj,knext,yield);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\n sigproj = " << SigProj << endl;
		str << "\n k = " << knext << endl;
		str << "\n yield = " << yield << endl; 
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	//// perform convergence verifications ---------------------------------------------------
	STATE theta = M_PI/5.;
	STATE xi = (Tensor2[0]+Tensor2[1]+Tensor2[2])/sqrt(3.);
	STATE beta = M_PI/3.;
	STATE k = -0.3;
	TPZManVector<STATE> xnorm,errnorm,converge;
	
	materialmodel.TaylorCheckDistF1(Tensor2, xi, beta, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck DistF1:" << endl;
    str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;		
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	materialmodel.TaylorCheckDDistF1DSigtrial(Tensor2, xi, beta, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck DDistF1DSigtrial:" << endl;
    str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//cout << "Teste de D(DistF2)/D(theta,beta,k)\n";
	{
		STATE kproj;
		TPZManVector<STATE,3> sigproj(3);
		materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
		STATE theta,beta;
		materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);
		materialmodel.TaylorCheckDistF2(Tensor2, theta, beta, kproj, kprev, xnorm, errnorm);
	}
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(DistF2)/D(theta,beta,k)" << endl;
    str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//cout << "Teste da derivada D(ResF2)/D(theta,beta,k)\n";
	materialmodel.TaylorCheckDDistF2(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(ResF2)/D(theta,beta,k):" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//cout << "Teste da derivada D(ResF2)/D(sigtrial)\n";
	{
		STATE kproj;
		TPZManVector<STATE,3> sigproj(3);
		materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
		STATE theta,beta;
		materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);
		materialmodel.TaylorCheckDDistF2DSigtrial(Tensor2, theta, beta, k, kprev, xnorm, errnorm);
	}
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(ResF2)/D(sigtrial):" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	materialmodel.TaylorCheckDF1Cart(xi, beta, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck DF1Cart:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//cout << "Teste de D(F2cart)/D(theta,beta,k)\n";
	{
		STATE kproj;
		TPZManVector<STATE,3> sigproj(3);
		materialmodel.ProjectF2(Tensor2, kprev, sigproj, kproj);
		STATE theta,beta;
		materialmodel.SurfaceParamF2(sigproj, kproj, theta, beta);
		materialmodel.TaylorCheckDF2Cart(theta, beta, kproj, xnorm, errnorm);
	}
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(F2cart)/D(theta,beta,k):" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//std::cout << "Teste de D(sigproj)/D(sigtrial) geral\n";
	materialmodel.TaylorCheckProjectSigma(Tensor2, kprev, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(sigproj)/D(sigtrial) geral:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	materialmodel.TaylorCheckParamF1Sigtrial(Tensor2,kprev,xnorm,errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck ParamF1Sigtrial:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	materialmodel.TaylorCheckProjectF1(Tensor2, kprev, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck ProjectF1:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//std::cout << "Teste de D(sigproj)/D(sigtrial) para funcao F2\n";
	materialmodel.TaylorCheckProjectF2(Tensor2, kprev, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(sigproj)/D(sigtrial) para funcao F2:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
	//std::cout << "Teste de D(theta,beta,k)/D(sigtrial) para funcao F2\n";
	materialmodel.TaylorCheckDtbkDsigtrial(Tensor2, kprev, xnorm, errnorm);
	materialmodel.ConvergenceRate(xnorm, errnorm, converge);
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\nTaylorCheck D(theta,beta,k)/D(sigtrial) para funcao F2:" << endl;
		str << "\n xnorm " << xnorm << endl;
    str << "\n errnorm " << errnorm << endl;
    str << "\n convergence rate " << converge << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
	
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
    TPZFMatrix<STATE> GradSigma(3,3);
    
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
        //materialmodel.ProjectSigma(sigtrial,kprev,sigproj,kproj);
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

void CheckDepConv()
{
	TPZSandlerExtended SDPV;
	SDPV.MCormicRanchSand(SDPV);
	STATE epsp=0.,k0;
	SDPV.Firstk(epsp, k0);
	TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV(k0);
	PlasticStepPV.fYC = SDPV;
	
	PlasticStepPV.fER = SDPV.GetElasticResponse();
	TPZElasticResponse ER = PlasticStepPV.fER; 
	STATE A=0.25,B=0.67,C=0.18,D=0.67,R=2.5,W=0.066,E=PlasticStepPV.fER.E(),nu=PlasticStepPV.fER.Poisson();
	
	TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick;
	PlasticStepErick.SetUp(nu, E, A, B, C, R,D,W);
	
	TPZTensor<STATE> eps1,eps2,eps3,sigma1,sigma2,sigma3,deps,stresstrial,sigmasol1,sigmasol2;
	
	TPZFNMatrix<36> DepPV(6,6,0.),DepER(6,6,0.);
	
	
	//    {0.0128547, -0.0564273, -0.0564273} beta=0
	//    {-0.0333333, 0.00666667, -0.0733333}beta=pi/2
//    STATE alpha = 1.e-3;
////	sigma1.XX()= 0.0128547;
//	sigma1.XX()=-0.0564273;
//	sigma1.YY()=-0.0564273;
//	sigma1.ZZ()=-0.0564273;
//    //sigma1.XX()=alpha*2;
//	//sigma1.YY()=-alpha*1;
//	//sigma1.ZZ()=-alpha*1;
//	//sigma1.XX()=-0.0333333;
//	//sigma1.YY()=0.00666667;
//	//sigma1.ZZ()= -0.0733333;
//	
//	ER.ComputeDeformation(sigma1, eps1);
	
/*	//    {0.121722, 0.0091389, 0.0091389}
	sigma1.XX()= 0.121722;
	sigma1.YY()= 0.0091389;
	sigma1.ZZ()= 0.0091389;
	
	ER.ComputeDeformation(sigma1, eps1);
	
	//    {0.135295, 0.0573526, 0.0573526}
	sigma1.XX()=0.135295;
	sigma1.YY()=0.0573526;
	sigma1.ZZ()=0.0573526;
	
	ER.ComputeDeformation(sigma1, eps1);
	
	TPZManVector<STATE> sigtr(3),sigproj(3);
	sigtr[0]=sigma1.XX();
	sigtr[1]=sigma1.YY();
	sigtr[2]=sigma1.ZZ();
	//STATE k1;
	//SDPV.ProjectF2(sigtr, k0,sigproj, k1);

 */
 
	//PlasticStepErick.ApplyStrainComputeDep(eps1,sigmasol1,DepER);
	//cout << " \n sigmasol1 = "<< sigmasol1 << endl;
	
	PlasticStepPV.ApplyStrainComputeDep(eps1, sigmasol2,DepPV);
	cout << " \n sigmasol2 "<< sigmasol2 << endl;
	
#ifdef LOG4CXX
	{
		std::stringstream outfile;
		outfile<<"\n Comparing Dep Matrix "<<std::endl;
		outfile<<"\n DepER "<< DepER << std::endl;
		outfile<<"\n DepPV "<< DepPV << std::endl;
		LOGPZ_DEBUG(logger,outfile.str());
		
	}
#endif
	
    
    
    sigma1.XX()= 0.0128547;
	//sigma1.XX()=-0.0564273;
	sigma1.YY()=-0.0564273;
	sigma1.ZZ()=-0.0564273;
    //sigma1.XX()=alpha*2;
	//sigma1.YY()=-alpha*1;
	//sigma1.ZZ()=-alpha*1;
	//sigma1.XX()=-0.0333333;
	//sigma1.YY()=0.00666667;
	//sigma1.ZZ()= -0.0733333;
	
	ER.ComputeDeformation(sigma1, eps1);
    TPZTensor<STATE> dsig;
    STATE alpha = 1.e-5;
    dsig.XX()= alpha;
    dsig.YY()=alpha;
    dsig.ZZ()=alpha;

    //dsig.XX()= 0;
    //dsig.YY()=alpha;
    //dsig.ZZ()=-alpha;
    
    //ER.ComputeDeformation(dsig, deps);
    
    //sigma1+=dsig;
	
#ifdef LOG4CXX
	{
		std::stringstream str;
		str << "\n-------------------------- BEGINING OF TAYLORCHECK TEST FOR DSIGDEPS USING PRINCIPAL VALUES --------------------------" << endl;

		TPZManVector<REAL,3> conv;
        deps.XX()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nXX:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
				
		deps.XX()=0.;
		deps.YY()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nYY:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		deps.YY()=0.;
		deps.ZZ()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nZZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		deps.ZZ()=0.;
		deps.XY()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nXY:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		deps.XY()=0.;
		deps.XZ()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nXZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;

		deps.XZ()=0.;
		deps.YZ()=-alpha;
		PlasticStepPV.TaylorCheck(eps1, deps, k0, conv);
		str << "\nYZ:" << endl;
		deps.Print(str);
		str << "conv = " << conv << endl;
		
		str << "\n-------------------------- END OF TAYLORCHECK TEST FOR DSIGDEPS USING PRINCIPAL VALUES --------------------------" << endl;
		LOGPZ_DEBUG(logger,str.str())
	}
#endif
	
}

void UniaxialLoadingPV(TPZPlasticStepPV<TPZYCMohrCoulombPV,TPZElasticResponse> PlasticStepPV)
{
	TPZTensor<STATE> epst,sigma1,sigload;
	TPZFNMatrix<36> DepPV(6,6,0.);
	STATE deltaeps = -0.001;
	epst.XX()=deltaeps;
	TPZElasticResponse ER;
	ER.SetUp(100, 0.25);
	ER.Compute(epst,sigload);
	//PlasticStepPV.ApplyStrainComputeDep(epst, sigma1,DepPV);
	//epst.XX()=0.;
	//PlasticStepPV.SetTensionSign(1);
	epst.XX()=0;
	sigload*=1;
	PlasticStepPV.ApplyLoad(sigload, epst);
///////	epst*=PlasticStepPV.SignCorrection();
	// PlasticStepPV.ApplyStrainComputeSigma(epst,sigma1);
}


void UniaxialLoadingPV(TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> PlasticStepPV,TPZVec<TPZPlasticState<REAL> > &state)
{
    
    

    TPZTensor<STATE> epst,sigma1;
    TPZFNMatrix<36> DepPV(6,6,0.);
    STATE deltaeps = -0.0005;
	epst.XX()=deltaeps;
    epst.YY()=deltaeps;
    epst.ZZ()=deltaeps;
    int nsteps=20;
    state.Resize(nsteps);
    for(int i=0;i<20;i++)
    {
        PlasticStepPV.ApplyStrainComputeDep(epst,sigma1,DepPV);
        state[i]=PlasticStepPV.GetState();
        epst.XX()+=deltaeps;
    }
    
/*//    //epst.YY()=deltaeps;
//    //epst.ZZ()=deltaeps;
//    TPZElasticResponse ER;
//    ER.SetUp(100, 0.25);
//    ER.Compute(epst,sigload);
//    PlasticStepPV.ApplyLoad(sigload, epst);
//    PlasticStepPV.ApplyStrainComputeSigma(epst, sig);
    dsig.XX()=-0.4/100.;
    dsig.YY()=dsig.XX()*0.4;
    dsig.ZZ()=dsig.XX()*0.4;
    ofstream outfile("Proportional04.txt");
    for(int i=0;i<100;i++)
    {

        PlasticStepPV.ApplyLoad(sigload, epst);
        outfile <<-epst.XX() << " " <<-sigload.XX()  << endl;
        cout << "\n state PV "<<PlasticStepPV.fN << endl;
        PlasticStepErick.ApplyLoad(sigload, epst2);
        cout << "\n state Erick "<<PlasticStepErick.GetState() << endl;
        
        sigload+=dsig;
        
    }
    
    PlasticStepErick.SetState(prevstate);
    PlasticStepPV.SetState(prevstate);
    
    
    dsig.XX()=-0.8/100.;
    dsig.YY()=dsig.XX()*0.6;
    dsig.ZZ()=dsig.XX()*0.6;
    for(int i=0;i<6;i++)sigload.fData[i]=0.;
    ofstream outfile2("Proportional06.txt");
    for(int i=0;i<100;i++)
    {
        PlasticStepPV.ApplyLoad(sigload, epst);
        outfile2 <<-epst.XX() << " " <<-sigload.XX()  << endl;
        cout << "\n state PV "<<PlasticStepPV.fN << endl;
        PlasticStepErick.ApplyLoad(sigload, epst2);
        cout << "\n state Erick "<<PlasticStepErick.GetState() << endl;
        sigload+=dsig;
        
    }
    
    
    PlasticStepErick.SetState(prevstate);
    PlasticStepPV.SetState(prevstate);
    
    dsig.XX()=-1.2/100.;
    dsig.YY()=dsig.XX()*0.8;
    dsig.ZZ()=dsig.XX()*0.8;
    for(int i=0;i<6;i++)sigload.fData[i]=0.;
    ofstream outfile3("Proportional08.txt");
    for(int i=0;i<100;i++)
    {
        PlasticStepPV.ApplyLoad(sigload, epst);
        outfile3 <<-epst.XX() << " " <<-sigload.XX()  << endl;
        cout << "\n state PV "<<PlasticStepPV.fN << endl;
        PlasticStepErick.ApplyLoad(sigload, epst2);
        cout << "\n state Erick "<<PlasticStepErick.GetState() << endl;
        sigload+=dsig;
        
    }
 */

}


void UniaxialLoadingErick(TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> PlasticStepErick,TPZVec<TPZPlasticState<REAL> > &state)
{
	
	TPZTensor<STATE> epst,sigma1;
	TPZFNMatrix<36> DepER(6,6,0.);
	STATE deltaeps = -0.0005;
	epst.XX()=deltaeps;
    epst.YY()=deltaeps;
    epst.ZZ()=deltaeps;
    int nsteps=20;
    state.Resize(nsteps);
    for(int i=0;i<nsteps;i++)
    {
        PlasticStepErick.ApplyStrainComputeDep(epst,sigma1,DepER);
        state[i]=PlasticStepErick.GetState();
        epst.XX()+=deltaeps;
        
    }
    
	
}

//template class TPZPlasticStepPV<TPZSandlerExtended, TPZElasticResponse>;


