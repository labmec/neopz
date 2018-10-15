//
//  main.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/20/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>
using namespace std;

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzinterpolationspace.h"

#include "pzgeoelbc.h"
#include "TPZCompElDisc.h"
#include "pzpoisson3d.h"
#include "pzelasmat.h"
#include "pzporoelasticmf2d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"

#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"

#include "mymeshes.h"
#include "pzfunction.h"

#include "pzlog.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <math.h>
using namespace std;


REAL const Pi = 4.*atan(1.);
REAL ftimeatual = 0.;
REAL fdeltaT = 0.;
REAL felarea = 0.;
REAL fbeta = 0.;
bool fdimensionless = false;


void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
void SolucaoPQMurad(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);//para calcular o erro
void SolucaoUMurad(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv);//para calcular o erro

void TerzaghiProblem1D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
void SolucaoPQTerzaghi(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

void BarryMercerProblem2D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
void SolucaoPQBarryMercer(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
void SolucaoUBarryMercer(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv);
void BarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);
void SolUBarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv);
void SolPBarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux);

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref);
void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area);

void ForcingSource(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingBCDeslocamento(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);


void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out,void (*fp)(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv));
void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out,void (*fp)(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv));


#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.porolasticmf2d.data"));
#endif

//problema murad e Loula
int main(int argc, char *argv[])
{
    
#ifdef LOG4CXX
	std::string logs("../logporoelastc2d.cfg");
	InitializePZLOG("../logporoelastc2d.cfg");
#endif
    
    bool triang=false;
    bool dimensionless = true;
    
    int nthreads = 0;
    
    REAL Eyoung = 3.e4;
    REAL poisson = 0.2;
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx=0.0;
    REAL fy = gravity*rockrho;
    
    REAL alpha = 1.0;
    REAL Se = 0.0;
    REAL perm = 1.e-10;
    REAL visc = 1.e-3;
    REAL sig0 = -1000.;
    REAL pini = 1000.;
    REAL La = 1.;
    REAL Ly = 5.*La;
    REAL Lx = 8.*La;
    REAL timeT = 10.;
    
    
    //Incompressible fluid (Se=0)
    if(dimensionless==true)
    {
        REAL lambda = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));//firstlame
        REAL mu = 0.5*Eyoung/(1+poisson);//secondlame
        REAL pref = pini;
        REAL Lref = La;
        REAL Se_aux = alpha*alpha/(lambda+2.*mu);
        REAL kovervisc = perm/visc;
        REAL Cf;//fluid diffusivity coeffcient
        REAL lambdaD = lambda*Se_aux;
        REAL muD = mu*Se_aux;
        Cf = kovervisc/Se_aux;

        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = sig0/pref;
        pini = pini/pref;
        timeT = timeT*Cf/(Lref*Lref);//1.
        Lx = Lx/Lref;
        Ly = Ly/Lref;
        perm = 1.;
        visc = 1.;
        fx = fx*Lref/pref;
        fy = fy*Lref/pref;
    }
    
    DadosMalhas * mydata = new DadosMalhas();
    mydata->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    
    ofstream saidaerro("ErroLoula.txt");
    
    
    for(int p = 2; p < 3; p++)
    {
        int pu = p;
        int pq = pu;
        int pp;
        if(triang==true){
            pp = pq-1;
        }else{
            pq=pu-1;
            pp = pq;
        }
        
        int h;
        saidaerro<<"\n CALCULO DO ERRO, ELEM. TRIANG., COM ORDEM POLINOMIAL pu = "<< pu << ", pq = "<< pq << " e pp = "<< pp<<endl;
        for (h = 3; h< 4; h++)
        {
        
        saidaerro<<"\n========= PARA h = "<< h<<"  ============= "<<endl;
    
        //primeira malha
        // geometric mesh (initial)
        //TPZGeoMesh * gmesh = mydata->GMesh(triang,Lx,Ly);
            
        TPZGeoMesh * gmesh = mydata->GMesh2(Lx,Ly,La);
        mydata->RefiningNearLine(2, gmesh, 0);
        mydata->AjustarContorno(gmesh);
            
        {
            std::ofstream malhaOut("malhageometrica.vtk");
            TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
            
            std::ofstream out("gmesh.txt");
            gmesh->Print(out);
        }

            
        //TPZGeoMesh * gmesh = mydata->GMesh3(LxD,LyD);
//        ofstream arg1("gmesh_inicial.txt");
//        gmesh->Print(arg1);

        // mydata->UniformRefine(gmesh, 2);


        // First computational mesh
        TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu,false,true);
//        ofstream arg2("cmesh1_inicial.txt");
//        cmesh1->Print(arg2);

        // second computational mesh
        TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq,false,true);
//        ofstream arg3("cmesh2_inicial.txt");
//        cmesh2->Print(arg3);


        // Third computational mesh
        TPZCompMesh * cmesh3 = mydata->CMeshPressure(gmesh, pp,triang,false);
        ofstream arg4("cmesh3_inicial.txt");
        cmesh3->Print(arg4);


//        // Cleaning reference of the geometric mesh to cmesh1
//        gmesh->ResetReference();
//        cmesh1->LoadReferences();
//        TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,h,false);
//        cmesh1->AdjustBoundaryElements();
//        cmesh1->CleanUpUnconnectedNodes();
////        ofstream arg5("cmesh1_final.txt");
////        cmesh1->Print(arg5);
//
//
//        // Cleaning reference to cmesh2
//        gmesh->ResetReference();
//        cmesh2->LoadReferences();
//        TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,h,false);
//        cmesh2->AdjustBoundaryElements();
//        cmesh2->CleanUpUnconnectedNodes();
////        ofstream arg6("cmesh2_final.txt");
////        cmesh2->Print(arg6);
//
//        // Cleaning reference to cmesh3
//        gmesh->ResetReference();
//        cmesh3->LoadReferences();
//        TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,h,true);
//
//        cmesh3->AdjustBoundaryElements();
//        cmesh3->CleanUpUnconnectedNodes();
////        ofstream arg7("cmesh3_final.txt");
////        cmesh3->Print(arg7);


        //	Set initial conditions for pressure

        TPZAnalysis an3(cmesh3,false);
        mydata->SolveSist(an3, cmesh3);
        int nrs = an3.Solution().Rows();
        TPZVec<STATE> solini(nrs,pini);
        //	cmesh3->Solution() = solucao1;
        //    cmesh3->LoadSolution(solucao1);

        TPZCompMesh  * cmeshL2 = mydata->CMeshPressureL2(gmesh, pp, solini,triang);
        TPZAnalysis anL2(cmeshL2,false);
        mydata->SolveSist(anL2, cmeshL2);
        //anL2.Solution().Print("sol");
        
        an3.LoadSolution(anL2.Solution());
        //    an3.Solution().Print();

        //malha multifisica
        TPZVec<TPZCompMesh *> meshvec(3);
        meshvec[0] = cmesh1;
        meshvec[1] = cmesh2;
        meshvec[2] = cmesh3;
        TPZPoroElasticMF2d * mymaterial;
        TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolucaoExata1D, 5);
        TPZCompMesh * mphysics = mydata->MalhaCompMultphysics(gmesh,meshvec,mymaterial,solExata);
//        ofstream arg8("mphysic.txt");
//        mphysics->Print(arg8);

        int NDeltaT = 1000000;
        int intervsaidas = NDeltaT/20;
        REAL deltaT=timeT/NDeltaT; //second
        mymaterial->SetTimeStep(deltaT);
        REAL maxTime = timeT;

       // mydata->SolveSistTransient(deltaT, maxTime, mymaterial, meshvec, mphysics,intervsaidas,ftimeatual);


        //Saida dos erros
//        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
//        TPZVec<REAL> erros;
//
//        saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
//        TPZAnalysis an12(cmesh1);
//        an12.SetExact(*SolucaoUMurad);
//        an12.PostProcessError(erros, saidaerro);
//
//
//        saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
//        TPZAnalysis an22(cmesh2);
//        an22.SetExact(*SolucaoPQMurad);
//        an22.PostProcessError(erros, saidaerro);
//
//        saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
//        TPZAnalysis an32(cmesh3);
//        an32.SetExact(*SolucaoUMurad);
//        an32.PostProcessError(erros, saidaerro);

        
        TPZAnalysis an(mphysics);
        TPZFMatrix<STATE> Initialsolution = an.Solution();
        //Initialsolution.Print("solini");
        
        std::string outputfile;
        outputfile = "TransientSolution";
        
//        std::stringstream outputfiletemp;
//        outputfiletemp << outputfile << ".vtk";
//        std::string plotfile = outputfiletemp.str();
//        mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);
        
        //Criando matriz de massa (matM)
        TPZAutoPointer <TPZMatrix<STATE> > matM = mydata->MassMatrix(mymaterial, mphysics,nthreads);
        
        //Criando matriz de rigidez (matK) e vetor de carga
        TPZFMatrix<STATE> matK;
        TPZFMatrix<STATE> fvec;
        mydata->StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec,nthreads);
        
        int nrows;
        nrows = matM->Rows();
        TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
        TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
        TPZFMatrix<STATE> Lastsolution = Initialsolution;
        
        REAL TimeValue = 0.0;
        int cent = 1;
        TimeValue = cent*deltaT;
        while (TimeValue <= maxTime)
        {
            ftimeatual  = TimeValue;
            // This time solution i for Transient Analytic Solution
            mymaterial->SetTimeValue(TimeValue);
            matM->Multiply(Lastsolution,TotalRhstemp);
            
            TotalRhs = fvec + TotalRhstemp;
            an.Rhs() = TotalRhs;
            an.Solve();
            Lastsolution = an.Solution();
            
            if(cent%intervsaidas==0){
                saidaerro<<"\n========= PARA O PASSO n = "<< cent <<"  E TEMPO tn = "<< TimeValue <<" =========\n"<<endl;
                std::stringstream outputfiletemp;
                outputfiletemp << outputfile << ".vtk";
                std::string plotfile = outputfiletemp.str();
                mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);
            
            
                TPZVec<REAL> erros;
                
                saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
                TPZAnalysis an12(cmesh1);
                an12.SetExact(*SolucaoUMurad);
                bool store_errors = false;
                an12.PostProcessError(erros, store_errors, saidaerro);
                
                
                saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
                TPZAnalysis an22(cmesh2);
                an22.SetExact(*SolucaoPQMurad);
                an22.PostProcessError(erros, store_errors, saidaerro);
                
                saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
                TPZAnalysis an32(cmesh3);
                an32.SetExact(*SolucaoPQMurad);
                an32.PostProcessError(erros, store_errors, saidaerro);
            }
            
            
            cent++;
            TimeValue = cent*deltaT;
        }

        cmesh1->CleanUp();
        cmesh2->CleanUp();
        cmesh3->CleanUp();
        //mphysics->CleanUp();
        delete cmesh1;
        delete cmesh2;
        delete cmesh3;
        //delete mphysics;
        delete gmesh;
    }
    }
    return 0;
}

void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
    
    bool sol_dimensionless=true;
    
	//REAL x = ptx[0];
	REAL x = ptx[1];
	
	REAL pini = 1000.;
	REAL lamb = 8333.33;
	REAL mi = 12500.0;
	REAL visc =0.001;
	REAL perm =  1.e-10;
	REAL H=1.;
	REAL tp = ftimeatual;
	int in;
	REAL pD = 0.0, uD = 0.0, sigD = 0.0, VDy=0.0;
    REAL sumpD = 0.0, sumuD = 0.0, sumsigD=0.0, sumVDy=0.0;
    
    REAL M =0.;
	REAL PI = atan(1.)*4.;
	
    sol.Resize(5, 0.);//p, ux, uy, sigx, sigy
    sol[0]=sol[1]=sol[2]=sol[3]=sol[4]=0.;
    flux.Resize(2, 1);
    flux(0,0)=0.;
    flux(1,0)=0.;
	
	REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
	REAL xD = fabs(1.-x)/H;
	for (in =999; in >= 0; in--) {
		
		M = PI*(2.*in+1.)/2.;
		sumpD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
		sumuD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
		sumsigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
        
        sumVDy+= 2.*cos(M*xD)*exp(-1.*M*M*tD);
	}
	
    pD = sumpD;
    uD = (H/H - xD) - sumuD;
    sigD = -1. + sumsigD;
    VDy = sumVDy;
    
    if(sol_dimensionless==true){
        sol[0] = pD;
        sol[2] = (-1.)*uD;
        sol[4] = sigD;
        
        flux(1,0) = (1.)*VDy;
    }else{
    
        sol[0] = pD*pini;
        sol[2] = (-1.)*uD*(pini*H)/(lamb+2.*mi);
        sol[4] = (sigD)*pini;
        flux(1,0) = (1.)*VDy*pini*(perm/visc);
    }
}


void SolucaoUMurad(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){
    
    bool sol_dimensionless=true;
    
	//REAL x = ptx[0];
	REAL x = ptx[1];
	
	REAL pini = 1000.;
	REAL lamb = 8333.33;
	REAL mi = 12500.0;
	REAL H=1.;
	REAL tp = ftimeatual;
	int in;
	REAL uD = 0.0, sigD=0.;
    REAL sumuD = 0.0, sumsigD = 0.;
    
    REAL M =0.;
	REAL PI = atan(1.)*4.;
	
    sol.Resize(2, 0.);// ux, uy;
    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;

	
	REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
	REAL xD = fabs(1.-x)/H;
	for (in =999; in >= 0; in--) {
		
		M = PI*(2.*in+1.)/2.;
		sumuD += (2./(M*M))*cos(M*xD)*exp(-1.*M*M*tD);
		sumsigD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
	}
	
    uD = (H/H - xD) - sumuD;
    sigD = -1. + sumsigD;
    
    if(sol_dimensionless==true){
        sol[1] = (-1.)*uD;
        deriv(1,1) = sigD;
    }else{
        sol[1] = (-1.)*uD*(pini*H)/(lamb+2.*mi);
        deriv(1,1) = (sigD)*pini;
    }
}


void SolucaoPQMurad(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
    
    bool sol_dimensionless=true;
    
	//REAL x = ptx[0];
	REAL x = ptx[1];
	
	REAL pini = 1000.;
	REAL visc =0.001;
	REAL perm =  1.e-10;
	REAL H=1.;
	REAL tp = ftimeatual;
	int in;
	REAL pD = 0.0, VDy=0.0;
    REAL sumpD = 0.0, sumVDy=0.0;
    
    REAL M =0.;
	REAL PI = atan(1.)*4.;
	
    sol.Resize(1, 0.);//p, ux, uy, sigx, sigy
    flux.Redim(3, 1);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
	
	REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
	REAL xD = fabs(1.-x)/H;
	for (in =999; in >= 0; in--) {
		
		M = PI*(2.*in+1.)/2.;
		sumpD += (2./M)*sin(M*xD)*exp(-1.*M*M*tD);
        sumVDy+= 2.*cos(M*xD)*exp(-1.*M*M*tD);
	}
	
    pD = sumpD;
    VDy = sumVDy;
    
    if(sol_dimensionless==true){
        sol[0] = pD;
        
        flux(1,0) = (1.)*VDy;
    }else{
        
        sol[0] = pD*pini;        
        flux(1,0) = (1.)*VDy*pini*(perm/visc);
    }
}



//problema de Terzaghi com Se!=0
int main_Terzaghi(int argc, char *argv[]){
#ifdef LOG4CXX
	std::string logs("../logporoelastc2d.cfg");
	InitializePZLOG("../logporoelastc2d.cfg");
#endif
    
    int nthreads = 0;
    
    bool triang = true;
    fdimensionless = true;
    
    REAL Eyoung = 1.e5;
    REAL poisson = 0.2;
//    REAL firstlame = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
//    REAL secondlame = 0.5*Eyoung/(1.+poisson);
    REAL alpha = 1.0;
    REAL Se = 1.e-1;
    
//    REAL K = firstlame + (2./3.)*secondlame;
    REAL sig0 = -1000.;
    REAL F  = abs(sig0);
    REAL pini = 0.;//(alpha*F)/(Se*(Ku+(4./3.)*secondlame));
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx=0.0;
    REAL fy = gravity*rockrho;
    REAL perm = 1.e-10;
    REAL visc = 1.e-4;
   
    REAL Ly = 1.;
    REAL Lx = 1.;
    REAL timeT = 100000;
    
    
    //Incompressible fluid (Se=0)
    if(fdimensionless==true)
    {
        REAL lambda = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));//firstlame
        REAL mu = 0.5*Eyoung/(1+poisson);//secondlame
        REAL lambdau = alpha*alpha/Se +lambda; //undrained;
        REAL Se_aux = alpha*alpha*(lambdau+2.*mu)/((lambdau - lambda)*(lambda+2.*mu));
        REAL lambdaD = lambda*Se_aux;
        REAL muD = mu*Se_aux;
        
        REAL pref = (alpha*F)/(Se*(lambdau + 2.*mu));
        REAL Lref = Ly;
        REAL kovervisc = perm/visc;
        REAL Cf;//fluid diffusivity coeffcient
        Cf = kovervisc/Se_aux;
        
        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = sig0/pref;
        pini =pini/pref;
        timeT =/* 0.1;//*/ timeT*Cf/(Ly*Ly);
        Lx = Lx/Lref;
        Ly = Ly/Lref;
        perm = 1.;
        visc = 1.;
        fx = fx*Lref/pref;
        fy = fy*Lref/pref;
        Se = Se/Se_aux;
    }
    
    DadosMalhas * mydata = new DadosMalhas();
    mydata->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    
    
    ofstream saidaerro("ErroTerzaghi.txt");
    
    
    for(int p =1; p<3; p++)
    {
        int pu = p;
        int pq = pu;
        int pp;
        if(triang==true){
            pp = pq-1;
        }else{
            pq=pu-1;
            pp = pq;
        }
        
        int h;
        saidaerro<<"\n CALCULO DO ERRO, ELEM. QUAD., COM ORDEM POLINOMIAL pu = "<< pu << ", pq = "<< pq << " e pp = "<< pp<<endl;
        for (h = 0; h< 6; h++)
        {
            
            saidaerro<<"\n========= PARA h = "<< h<<"  ============= "<<endl;
            // geometric mesh (initial)
            TPZGeoMesh * gmesh = mydata->GMesh(triang,Lx,Ly);
            mydata->UniformRefine(gmesh, h);
            //mydata->RefiningNearLine(2, gmesh, 4);


            // First computational mesh
            TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu,false,false);

            // second computational mesh
            TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq, false,false);

            // Third computational mesh
            TPZCompMesh * cmesh3=mydata->CMeshPressure(gmesh, pp,triang, false);
/*
            // Cleaning reference of the geometric mesh to cmesh1
            gmesh->ResetReference();
            cmesh1->LoadReferences();
            TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,h,false);
            cmesh1->AdjustBoundaryElements();
            cmesh1->CleanUpUnconnectedNodes();


            // Cleaning reference to cmesh2
            gmesh->ResetReference();
            cmesh2->LoadReferences();
            TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,h,false);
            cmesh2->AdjustBoundaryElements();
            cmesh2->CleanUpUnconnectedNodes();

            // Cleaning reference to cmesh3
            gmesh->ResetReference();
            cmesh3->LoadReferences();
            TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,h,true);

            cmesh3->AdjustBoundaryElements();
            cmesh3->CleanUpUnconnectedNodes();
*/

            //	Set initial conditions for pressure
            TPZAnalysis an3(cmesh3);
            int nrs = an3.Solution().Rows();
            TPZVec<STATE> solini(nrs,pini);
            TPZCompMesh  * cmeshL2 = mydata->CMeshPressureL2(gmesh, pp, solini,triang);
            TPZAnalysis anL2(cmeshL2);
            mydata->SolveSist(anL2, cmeshL2);
            //anL2.Solution().Print("sol");
            an3.LoadSolution(anL2.Solution());

            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(3);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            meshvec[2] = cmesh3;
            TPZPoroElasticMF2d * mymaterial;
            TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(TerzaghiProblem1D, 5);
            TPZCompMesh * mphysics = mydata->MalhaCompTerzaghi(gmesh,meshvec,mymaterial,solExata);

            
            int NDeltaT = 1000000;
            int intervsaidas = NDeltaT/20;
            REAL deltaT=timeT/NDeltaT; //second
            mymaterial->SetTimeStep(deltaT);
//            REAL maxTime = timeT;

            //mydata->SolveSistTransient(deltaT, maxTime, mymaterial, meshvec, mphysics,100,ftimeatual);

            
            //==== Imprimir erros ======
            TPZAnalysis an(mphysics);
            TPZFMatrix<STATE> Initialsolution = an.Solution();

            std::string outputfile;
            outputfile = "TransientSolution";

//        std::stringstream outputfiletemp;
//        outputfiletemp << outputfile << ".vtk";
//        std::string plotfile = outputfiletemp.str();
//        mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);

            //Criando matriz de massa (matM)
            TPZAutoPointer <TPZMatrix<STATE> > matM = mydata->MassMatrix(mymaterial, mphysics, nthreads);

            //Criando matriz de rigidez (matK) e vetor de carga
            TPZFMatrix<STATE> matK;
            TPZFMatrix<STATE> fvec;
            mydata->StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec, nthreads);

            int nrows;
            nrows = matM->Rows();
            TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
            TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
            TPZFMatrix<STATE> Lastsolution = Initialsolution;

            REAL TimeValue = 0.0;
            int cent = 1;
            TimeValue = cent*deltaT;
            //while (TimeValue <= maxTime)
            while (cent <= intervsaidas)
            {
                ftimeatual  = TimeValue;
                // This time solution i for Transient Analytic Solution
                mymaterial->SetTimeValue(TimeValue);
                matM->Multiply(Lastsolution,TotalRhstemp);
                
                TotalRhs = fvec + TotalRhstemp;
                an.Rhs() = TotalRhs;
                an.Solve();
                Lastsolution = an.Solution();
                
                if(cent%intervsaidas==0){
                    saidaerro<<"\n========= PARA O PASSO n = "<< cent <<"  E TEMPO tn = "<< TimeValue <<" =========\n"<<endl;
                    
//                    std::stringstream outputfiletemp;
//                    outputfiletemp << outputfile << ".vtk";
//                    std::string plotfile = outputfiletemp.str();
//                    mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);
                    
                    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
                    TPZVec<REAL> erros;
                    
                    saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
//                    TPZAnalysis an22(cmesh2);
//                    an22.SetExact(*SolucaoPQTerzaghi);
//                    an22.PostProcessError(erros, saidaerro);
                    
                    ErrorHDiv2(cmesh2,saidaerro,SolucaoPQTerzaghi);
                    
                    saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
//                    TPZAnalysis an32(cmesh3);
//                    an32.SetExact(*SolucaoPQTerzaghi);
//                    an32.PostProcessError(erros, saidaerro);
                    ErrorH1(cmesh3,saidaerro,SolucaoPQTerzaghi);
                }
                
                
                cent++;
                TimeValue = cent*deltaT;
            }

            cmesh1->CleanUp();
            cmesh2->CleanUp();
            cmesh3->CleanUp();
            //mphysics->CleanUp();
            delete cmesh1;
            delete cmesh2;
            delete cmesh3;
            //delete mphysics;
            delete gmesh;

        }
    }
    
    return 0;
}

void TerzaghiProblem1D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
    
    //REAL x = ptx[0];
	REAL y = 1.- ptx[1];
	
    REAL alpha = 1.;
    REAL Se = 1.e-1;
    
    REAL Eyoung = 1.e5;
    REAL poisson = 0.2;
    REAL lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
    REAL mu = 0.5*Eyoung/(1.+poisson);
    
	REAL visc =1.e-4;
	REAL perm =  1.e-10;
    
    REAL H=1.;
    REAL F = 1000.;
	REAL tp = ftimeatual;
	
    REAL kappa = perm/visc;
    REAL K = lamb + (2./3.)*mu;
    REAL Ku = K + alpha*alpha/Se;
    REAL Cf = (1./Se)*kappa*(K + (4./3.)*mu)/(Ku + (4./3.)*mu);
    REAL S0 = (F*H)/(Ku + (4./3.)*mu);
    REAL SInfty = (F*H)/(K + (4./3.)*mu);
    
    
	int in;
	REAL p = 0.0, uyH = 0.0, VDy = 0.0;
    REAL sump = 0.0, sumuyH = 0.0, sumVDy=0.0;
    
	REAL PI = atan(1.)*4.;
	
	sol.Resize(5, 0.);
    flux(0,0)=0.;
    flux(1,0)=0.;
	
    REAL aux1= (alpha*F)/(Se*(Ku + (4./3.)*mu));
    REAL aux2=0.;
    REAL aux3=0.;
    
    if(fdimensionless==true)
    {
        for (in =999; in >= 0; in--)
        {
            aux2 = 4./(PI*(2.*in + 1.));
            aux3 = PI*(2.*in + 1.)/2.;
            
            sump += aux2*sin(aux3*y)*exp(-aux3*aux3*tp);
            
            sumuyH += 0.5*aux2*aux2*exp(-aux3*aux3*tp);
            
            sumVDy += cos(aux3*y)*exp(-aux3*aux3*tp);
        }
        
        p = sump;
        uyH = SInfty + (S0 - SInfty)*sumuyH/aux1;
        VDy = (-2.)*sumVDy;
        
        sol[0] = p;
        sol[2] = -uyH;
        flux(1,0) = -VDy;
        
    }else{
        for (in =999; in >= 0; in--)
        {
            aux2 = 4./(PI*(2.*in + 1.));
            aux3 = Pi*(2.*in + 1.)/(2.*H);
            
            sump += aux2*sin(aux3*y)*exp(-aux3*aux3*Cf*tp);
            
            sumuyH += 0.5*aux2*aux2*exp(-aux3*aux3*Cf*tp);
            
            sumVDy += cos(aux3*y)*exp(-aux3*aux3*Cf*tp);
        }
        
        p = aux1*sump;
        uyH = SInfty + (S0 - SInfty)*sumuyH;
        VDy = -(2./H)*kappa*aux1*sumVDy;
        
        
        sol[0] = p;
        sol[2] = -uyH;
        flux(1,0) = -VDy;
    }
}

void SolucaoPQTerzaghi(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux)
{
 
    bool dimensionless = true;
    
    //REAL x = ptx[0];
	REAL y = 1.- ptx[1];
	
    REAL alpha = 1.;
    REAL Se = 1.e-1;
    
    REAL Eyoung = 1.e5;
    REAL poisson = 0.2;
    REAL lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
    REAL mu = 0.5*Eyoung/(1.+poisson);
    
	REAL visc =1.e-4;
	REAL perm =  1.e-10;
    
    REAL H=1.;
    REAL F = 1000.;
	REAL tp = ftimeatual;
	
    REAL kappa = perm/visc;
    REAL K = lamb + (2./3.)*mu;
    REAL Ku = K + alpha*alpha/Se;
    REAL Cf = (1./Se)*kappa*(K + (4./3.)*mu)/(Ku + (4./3.)*mu);
    
	int in;
	REAL p = 0.0, VDy = 0.0;// DivVD = 0.0;
    REAL sump = 0.0, sumVDy=0.0;//, sumDiv = 0.0;
    
	REAL PI = atan(1.)*4.;
	
	sol.Resize(1, 0.);
    flux.Resize(2,1);
    flux(0,0)=0.;
    flux(1,0)=0.;
    //flux(2,0)=0.;
	
    REAL aux1= (alpha*F)/(Se*(Ku + (4./3.)*mu));
    REAL aux2=0.;
    REAL aux3=0.;
    
    if(dimensionless==true)
    {
        for (in =999; in >= 0; in--)
        {
            aux2 = 4./(PI*(2.*in + 1.));
            aux3 = PI*(2.*in + 1.)/2.;
            
            sump += aux2*sin(aux3*y)*exp(-aux3*aux3*tp);
            
            sumVDy += cos(aux3*y)*exp(-aux3*aux3*tp);
            
            //sumDiv += (-aux3)*sin(aux3*y)*exp(-aux3*aux3*tp);
        }
        
        p = sump;
        VDy = -2.*sumVDy;
        //DivVD = 2.*sumDiv;
        
        sol[0] = p;
        flux(1,0) = -VDy;
        //flux(2,0) = -DivVD;
        
    }else{
        for (in = 999; in >= 0; in--)
        {
            aux2 = 4./(PI*(2.*in + 1.));
            aux3 = Pi*(2.*in + 1.)/(2.*H);
            
            sump += aux2*sin(aux3*y)*exp(-aux3*aux3*Cf*tp);
            
            sumVDy += cos(aux3*y)*exp(-aux3*aux3*Cf*tp);
            
            //sumDiv += (-aux3)*sin(aux3*y)*exp(-aux3*aux3*tp);
        }
        
        p = aux1*sump;
        VDy = (2./H)*kappa*aux1*sumVDy;
       // DivVD = (2./H)*kappa*aux1*sumDiv;
        
        sol[0] = p;
        flux(1,0) = -VDy;
        //flux(2,0) = -DivVD;
    }
}

//Problema Barry and Mercer
int main_BarryMercer(int argc, char *argv[]){

    #ifdef LOG4CXX
        InitializePZLOG();
    #endif
    
    gRefDBase.InitializeUniformRefPattern(EOned);
    gRefDBase.InitializeUniformRefPattern(EQuadrilateral);
    
    int nthreads = 0;
    
    bool triang = false;
    fdimensionless = true;
    
    REAL Eyoung = 1.e5;
    REAL poisson = 0.1;
    
    REAL alpha = 1.0;
    REAL Se = 0.;
    
    REAL sig0 =0.;
    REAL pini = 0.;
    
    REAL Lx = 1.;
    REAL Ly = 1.;
    
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx=0.0;
    REAL fy = gravity*rockrho;
    REAL perm = 1.e-6;
    REAL visc = 1.e-4;
    
    //trabalho com o tempo adimensional
    REAL timeT = 0.00153589;//0.00614356;
    REAL lambda = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));//firstlame
    REAL mu = 0.5*Eyoung/(1+poisson);//secondlame
    REAL Se_aux = alpha*alpha/(lambda+2.*mu);
    REAL kovervisc = perm/visc;
    REAL Cf = kovervisc/Se_aux;//fluid diffusivity coeffcient
    REAL beta = Cf/(Ly*Ly);
    fbeta = beta;
    timeT = timeT*beta;
    
    if(fdimensionless==true)
    {
        REAL Lref = sqrt(Lx*Ly);
        REAL lambdaD = lambda*Se_aux;
        REAL muD = mu*Se_aux;
        
        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = 0.;
        pini = 0.;
        Lx = Lx/Lref;
        Ly = Ly/Lref;
        perm = 1.;
        visc = 1.;
        fx = 0.;
        fy = 0.;
        Se = 0.;//Incompressible fluid (Se=0)
    }
    
    {//NAO APAGAR!!!
        char buf[] =
        "16 10 "
        "-50 Qua000022224 "
        "-1 -1 0 "
        "1 -1 0 "
        "1 1 0  "
        "-1 1 0 "
        "-0.33333 -1 0  "
        "+0.33333 -1 0  "
        "+1 -0.33333 0  "
        "+1 +0.33333 0  "
        "+0.33333 +1 0  "
        "-0.33333 +1 0  "
        "-1 +0.33333 0  "
        "-1 -0.33333 0  "
        "-0.33333 -0.33333 0    "
        "+0.33333 -0.33333 0    "
        "+0.33333 +0.33333 0    "
        "-0.33333 +0.33333 0    "
        "3 4 0  1  2  3 "
        "3 4 0 4 12 11  "
        "3 4 4 5 13 12  "
        "3 4 5 1 6 13   "
        "3 4 11 12 15 10    "
        "3 4 12 13 14 15    "
        "3 4 13 6 7 14  "
        "3 4 10 15 9 3  "
        "3 4 15 14 8 9  "
        "3 4 14 7 2 8   ";
        std::istringstream str(buf);
        TPZAutoPointer<TPZRefPattern> refp = new TPZRefPattern(str);
        refp->GenerateSideRefPatterns();
        gRefDBase.InsertRefPattern(refp);
        if(!refp)
        {
            DebugStop();
        }
    }
    
//    {//NAO APAGAR!!!
//        char buf[] =
//        "4 4 "
//        "-50 Lin000022224 "
//        "-1 0 0 "
//        "1 0 0 "
//        "-0.33333 0 0  "
//        "+0.33333 0 0  "
//        "1 2 0 1  "
//        "1 2 0 2  "
//        "1 2 2 3  "
//        "1 2 3 1  ";
//        std::istringstream str(buf);
//        TPZAutoPointer<TPZRefPattern> refp = new TPZRefPattern(str);
//        gRefDBase.InsertRefPattern(refp);
//        if(!refp)
//        {
//            DebugStop();
//        }
//    }

    
    DadosMalhas * mydata = new DadosMalhas();
    
    TPZVec<REAL> pt(3);
    pt[0] = Lx/4.;
    pt[1] = Ly/4.;
    pt[2] = 0.;
    REAL Area;
    
    ofstream saidaerro("Erro.txt");
    for(int p = 2; p<3; p++)
    {
        int pu = p+1;
        int pq = p;
        int pp = p;
        int h;
        
        saidaerro<<"\n CALCULO DO ERRO, ELEM. TRIANG., COM ORDEM POLINOMIAL pu = "<< pu << ", pq = "<< pq << " e pp = "<< pp<<endl;
        for (h = 1; h< 2; h++)
        {
            
            saidaerro<<"\n========= PARA h = "<< h<<"  ============= "<<endl;

            TPZGeoMesh * gmesh = mydata->GMesh4(Lx,Ly,0,0);
            mydata->UniformRefine(gmesh, h);
            RefinamentoPadrao3x3(gmesh,h+1);
            
            int newmatId = mydata->GetIdSourceTerm();//mat id of the source term
            RefinamentoPadrao3x3(gmesh,h+1,pt, true, newmatId, Area);
            felarea = Area;
            
//            std::ofstream malhaGeo("gmesh2D.txt");
//            gmesh->Print(malhaGeo);
            
            
            mydata->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
            
            
            // First computational mesh
            TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu,true,false);
            ofstream arg1("cmesh1.txt");
            cmesh1->Print(arg1);
            
            // second computational mesh
            TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq,true,false);
            ofstream arg2("cmesh2.txt");
            cmesh2->Print(arg2);
            
            // Third computational mesh
            TPZCompMesh * cmesh3=mydata->CMeshPressure(gmesh, pp,triang,true);
            ofstream arg3("cmesh3.txt");
            cmesh3->Print(arg3);

            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(3);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            meshvec[2] = cmesh3;
            
            TPZAutoPointer<TPZFunction<STATE> > sourceterm = new TPZDummyFunction<STATE>(ForcingSource, 5);
            TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(BarryMercerProblem2D, 5);
            
            TPZCompMesh * mphysics = mydata->MalhaCompBarryMercer(gmesh,meshvec,sourceterm, solExata);
            ofstream arg8("mphysic.txt");
            mphysics->Print(arg8);
            
//            ofstream arg9("gmesh_final.txt");
//            gmesh->Print(arg9);
            
            
            int NDeltaT = 100;
            int intervsaidas = NDeltaT/10;
            REAL deltaT=timeT/NDeltaT; //second
//            REAL maxTime = timeT;
            fdeltaT = deltaT;
            
            //mydata->SolveSistBarryMercert(deltaT, maxTime, meshvec, mphysics,1,ftimeatual);
            
            
            ///======== Impimir Erros ========
            TPZMaterial *mat1 = mphysics->FindMaterial(mydata->GetMatId());
            TPZMaterial *mat2 = mphysics->FindMaterial(mydata->GetMatId()+1);
            TPZPoroElasticMF2d * mat12 = dynamic_cast<TPZPoroElasticMF2d *>(mat1);
            TPZPoroElasticMF2d * mat22 = dynamic_cast<TPZPoroElasticMF2d *>(mat2);
            
            mat12->SetTimeStep(deltaT);
            mat22->SetTimeStep(deltaT);
            
            TPZAnalysis an(mphysics);
            TPZFMatrix<STATE> Initialsolution = an.Solution();
            
            std::string outputfile;
            outputfile = "TransientSolution";
            
            //Criando matriz de massa (matM)
            TPZAutoPointer <TPZMatrix<STATE> > matM = mydata->MassMatrix(mphysics, nthreads);
            
            //Criando matriz de rigidez (matK) e vetor de carga
            TPZFMatrix<STATE> matK;
            TPZFMatrix<STATE> fvec;
            mydata->StiffMatrixLoadVec(mphysics, an, matK, fvec, nthreads);
            
            int nrows;
            nrows = matM->Rows();
            TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
            TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
            TPZFMatrix<STATE> Lastsolution = Initialsolution;
            
            REAL TimeValue = 0.0;
            int cent = 1;
            TimeValue = cent*deltaT;
            while (cent<NDeltaT+1)//TimeValue < maxTime || TimeValue == maxTime)
            {
                ftimeatual  = TimeValue;
                // This time solution i for Transient Analytic Solution
                matM->Multiply(Lastsolution,TotalRhstemp);
                mydata->StiffMatrixLoadVec(mphysics, an, matK, fvec, nthreads);
                
                TotalRhs = fvec + TotalRhstemp;
                an.Rhs() = TotalRhs;
                an.Solve();
                Lastsolution = an.Solution();
                
                if(cent%intervsaidas==0)
                {
                    saidaerro<<"\n========= PARA O PASSO n = "<< cent <<"  E TEMPO tn = "<< TimeValue <<" =========\n"<<endl;
                    
                    std::stringstream outputfiletemp;
                    outputfiletemp << outputfile << ".vtk";
                    std::string plotfile = outputfiletemp.str();
                    mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);
//                    
//                    TPZVec<REAL> erros;
//            
//                    saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
//                    TPZAnalysis an12(cmesh1);
//                    an12.SetExact(*SolucaoUBarryMercer);
//                    an12.PostProcessError(erros, saidaerro);
//                    
//                    
//                    saidaerro<<" \nErro da simulacao multifisica do fluxo (q)" <<endl;
//                    TPZAnalysis an22(cmesh2);
//                    an22.SetExact(*SolucaoPQBarryMercer);
//                    an22.PostProcessError(erros, saidaerro);
//                    
//                    saidaerro<<" Erro da simulacao multifisica da pressao (p)" <<endl;
//                    TPZAnalysis an32(cmesh3);
//                    an32.SetExact(*SolucaoPQBarryMercer);
//                    an32.PostProcessError(erros, saidaerro);
                }
                
                cent++;
                TimeValue = cent*deltaT;
//                if(cent==200 || TimeValue == maxTime){
//                    TimeValue = cent*deltaT;
//                    int res = cent%intervsaidas;
//                }
            
            }
            
            cmesh1->CleanUp();
            cmesh2->CleanUp();
            cmesh3->CleanUp();
            //mphysics->CleanUp();
            delete cmesh1;
            delete cmesh2;
            delete cmesh3;
            //delete mphysics;
            delete gmesh;
        }
    }
    
    return 0;
}

void ForcingSource(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    REAL tp = ftimeatual;
//    REAL deltaT = fdeltaT;
    REAL elarea = felarea;
    
    //REAL temp =(cos(tp-deltaT) - cos(tp))/deltaT;
//    disp[0]= 2.*temp/deltaT;
    REAL temp = sin(tp);
    
    if(fdimensionless==true){
        disp[0]= 2.*temp/elarea;
    }
    else disp[0]= 2.*fbeta*temp/elarea;
}

void BarryMercerProblem2D(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
        
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL x0 = 0.25;
    REAL y0 = 0.25;
    REAL tp = ftimeatual;
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
   
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.1;
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;

    sol.Resize(5, 0.);
	
    flux(0,0)=0.;
    flux(1,0)=0.;
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil = 0., util = 0., wtil =0., temp=0.;
    REAL sump=0., sumux = 0., sumuy = 0., sumVDx = 0., sumVDy=0.;
    
    for(in=99; in<=0; in--)
    {
        gamman = in*PI;
        for(jq=99; jq<=0; jq--)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            temp = -2./(gammanq*gammanq+1.);
            
            ptil = temp*sin(gamman*x0)*sin(gammaq*y0)*(gammanq*sin(tp) - cos(tp) + exp(-gammanq*tp));
            util = gamman*ptil/gammanq;
            wtil = gammaq*ptil/gammanq;
            temp = 4./(Lx*Ly);
            
            sump += temp*ptil*sin(gamman*x)*sin(gammaq*y); //(multiplicado por -1, multipliquei depois da soma)
            sumux += temp*util*cos(gamman*x)*sin(gammaq*y);
            sumuy += temp*wtil*sin(gamman*x)*cos(gammaq*y);
            sumVDx += temp*ptil*gamman*cos(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
            sumVDy += temp*ptil*gammaq*sin(gamman*x)*cos(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
        }
    }
    
    if(fdimensionless == false){
        sol[0] = -(lamb+2.*mu)*sump;
        sol[1] = Lref*sumux;
        sol[2] = Lref*sumuy;
        flux(0,0) = kovervisc*(lamb+2.*mu)*sumVDx; //(multiplicado por -kovervisc, cancela com o outro sinal de menos)
        flux(1,0) = kovervisc*(lamb+2.*mu)*sumVDy; //(multiplicado por -kovervisc, cancela com o outro sinal de menos)
    }
    else{
        sol[0] = -sump;
        sol[1] = sumux;
        sol[2] = sumuy;
        flux(0,0) = sumVDx;//(multiplicado por -kovervisc, cancela com o outro sinal de menos)
        flux(1,0) = sumVDy;//(multiplicado por -kovervisc, cancela com o outro sinal de menos)
    }
}

void SolucaoUBarryMercer(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){
  
    
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL x0 = 0.25;
    REAL y0 = 0.25;
    REAL tp = ftimeatual;
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
    
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.1;
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(2, 0.);// ux, uy;
    deriv.Resize(2,2);//sigx, sigxy, sigyx, sigy
    deriv(0,0) = deriv(0,1) = deriv(1,0) = deriv(1,1) = 0.;
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil = 0., util = 0., wtil =0., temp=0.;
    REAL sumux = 0., sumuy = 0., sumsigx=0., sumsigy=0., sumsigxy=0.;
    
    for(in=499; in<=0; in--)
    {
        gamman = in*PI;
        for(jq=499; jq<=0; jq--)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            temp = -2./(gammanq*gammanq+1);
            
            ptil = temp*sin(gamman*x0)*sin(gammaq*y0)*(gammanq*sin(tp) - cos(tp) + exp(-gammanq*tp));
            util = gamman*ptil/gammanq;
            wtil = gammaq*ptil/gammanq;
            temp = 4./(Lx*Ly);
            
            sumux += temp*util*cos(gamman*x)*sin(gammaq*y);
            sumuy += temp*wtil*sin(gamman*x)*cos(gammaq*y);
            sumsigx += gamman*temp*util*sin(gamman*x)*sin(gammaq*y);//multiplicado por -1
            sumsigy += gammaq*temp*wtil*sin(gamman*x)*sin(gammaq*y);//multiplicado por -1
            sumsigxy += gammaq*temp*util*cos(gamman*x)*cos(gammaq*y);
        }
    }
    
    if(fdimensionless == false){
        std::cout<<"nao foi calculdo\n";
        DebugStop();
    }
    else{
        sol[0] = sumux;
        sol[1] = sumuy;
        deriv(0,0) = -sumsigx;
        deriv(1,1) = -sumsigy;
        deriv(0,1) = sumsigxy;
        deriv(1,0) = sumsigxy;
    }
    
}

void SolucaoPQBarryMercer(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
   
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL x0 = 0.25;
    REAL y0 = 0.25;
    REAL tp = ftimeatual;
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
    
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.1;
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(1, 0.);
	flux.Redim(3, 1);
    flux(0,0)=flux(1,0)=flux(2,0)=0.;
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil = 0., util = 0., wtil =0., temp=0.;
    REAL sump=0., sumVDx = 0., sumVDy=0.;
    
    for(in=499; in<=0; in--)
    {
        gamman = in*PI;
        for(jq=499; jq<=0; jq--)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            temp = -2./(gammanq*gammanq+1);
            
            ptil = temp*sin(gamman*x0)*sin(gammaq*y0)*(gammanq*sin(tp) - cos(tp) + exp(-gammanq*tp));
            util = gamman*ptil/gammanq;
            wtil = gammaq*ptil/gammanq;
            temp = 4./(Lx*Ly);
            
            sump += temp*ptil*sin(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
            sumVDx += temp*ptil*gamman*cos(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
            sumVDy += temp*ptil*gammaq*sin(gamman*x)*cos(gammaq*y);
        }
    }
    
    if(fdimensionless == false){
        sol[0] = -(lamb+2.*mu)*sump;
        flux(0,0) = kovervisc*(lamb+2.*mu)*sumVDx;//(multiplicado por -kovervisc, cancela com o outro sinal de menos)
        flux(1,0) = kovervisc*(lamb+2.*mu)*sumVDy;
    }
    else{
        sol[0] = -sump;
        flux(0,0) = sumVDx;//(multiplicado por -kovervisc, cancela com o outro sinal de menos)
        flux(1,0) = sumVDy;
    }

}


//void BarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
//    
//    REAL x = ptx[0];
//	REAL y = ptx[1];
//    
//    //dados do problema
//	REAL Lx = 1.;
//    REAL Ly = 1.;
//    REAL eps = 0.0;
//    REAL x0 = 0. + eps;
//    REAL x1 = Lx - eps;
//    REAL m = 1.5;
//    REAL tp = ftimeatual;
//    REAL Ft = sin(tp);
//    
//    REAL perm =0., visc=0., kovervisc=0.;
//    
//    //dados adimensional
//    if(fdimensionless == true)
//    {
//        perm = 1.;
//        visc = 1.;
//        kovervisc =perm/visc;
//    }
//    
//    
//    //dados com dimensao
//    REAL Lref = 0.;
//    REAL lamb=0., mu=0.;
//    if(fdimensionless == false){
//        
//        REAL Eyoung = 1.e5;
//        REAL poisson = 0.5*(1.-1./m);
//        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
//        mu = 0.5*Eyoung/(1.+poisson);
//        perm =  1.e-6;
//        visc =1.e-4;
//        kovervisc =perm/visc;
//        Lref = sqrt(Lx*Ly);
//        //beta = (lamb+2.*mu)*kovervisc/Lref;
//        //tp = tp;
//    }
//	
//    int in, jq;
//    REAL PI = atan(1.)*4.;
//    
//    sol.Resize(5, 0.);
//	
//    flux(0,0)=0.;
//    flux(1,0)=0.;
//    
//    REAL gamman = 0., gammaq = 0., gammanq = 0.;
//	REAL ptil = 0., util = 0., wtil =0.;
//    REAL sump=0., sumux = 0., sumuy = 0., sumVDx = 0., sumVDy=0., temp=0.;
//    
//    for(in=1; in<=100; in++)
//    {
//        gamman = in*PI;
//        temp = (eps*eps*eps)*(gamman*gamman*gamman*gamman);
//        
//        //Funcao Heaviside
////        ptil = ((cos(gamman*x0)-cos(gamman*x1)))/gamman;
//        
//        //Funcao Heaviside suave
////        if(x>=(x0 - eps) && x<=x0) ptil = (eps*gamman*(-(6.+ eps*eps*gamman*gamman)*cos(gamman*x0) - 6.*cos(gamman*(x0-eps))) + 12.*(sin(gamman*(eps - x0)) + sin(gamman*x0)))/temp;
////        
////        if(x>=x0 && x<=x1) ptil = ((cos(gamman*x0)-cos(gamman*x1)))/gamman;
////        
////        if(x>=x1 && x<=(x1+eps)) ptil = (eps*gamman*(6.+ eps*eps*gamman*gamman)*cos(gamman*x1) +
////            6.*eps*gamman*cos(gamman*(eps+x1)) + 12.*sin(gamman*x1) - 12.*sin(gamman*(eps+x1)))/temp;
////        
//         if(fabs(y)==0)
//         {
//             
//             //Funcao Heaviside suave
////            if(x>=(x0 - eps) && x<x0) ptil = (eps*gamman*(-(6.+ eps*eps*gamman*gamman)*cos(gamman*x0) - 6.*cos(gamman*(x0-eps))) + 12.*(sin(gamman*(eps - x0)) + sin(gamman*x0)))/temp;
////
////            if(x>=x0 && x<=x1) ptil = ((cos(gamman*x0)-cos(gamman*x1)))/gamman;
////
////            if(x>x1 && x<=(x1+eps)) ptil = (eps*gamman*(6.+ eps*eps*gamman*gamman)*cos(gamman*x1) +
////                6.*eps*gamman*cos(gamman*(eps+x1)) + 12.*sin(gamman*x1) - 12.*sin(gamman*(eps+x1)))/temp;
////            
////             ptil = ptil*Ft;
//             
//             ptil = ((cos(gamman*x0)-cos(gamman*x1))*Ft)/gamman;
//             sump += ptil*sin(gamman*x);
//             sumVDx += ptil*gamman*cos(gamman*x);//(multiplicado por -1, multipliquei depois da soma)
//             sumVDy = 0.;//(multiplicado por -1, multipliquei depois da soma)
//             
//             for(jq=1; jq<=100; jq++)
//             {
//                 gammaq = jq*PI;
//                 gammanq = gamman*gamman + gammaq*gammaq;
//                 
//                 util = (gamman*ptil)/gammanq;
//                 wtil = (gamman*gamman*ptil)/(gammanq*gammaq);
//                 
//                 sumux += util*cos(gamman*x)*sin(gammaq*y);
//                 sumuy += wtil*sin(gamman*x)*cos(gammaq*y);
//             }
//        
//         }else
//         {
//            for(jq=1; jq<=100; jq++)
//            {
//                gammaq = jq*PI;
//                gammanq = gamman*gamman + gammaq*gammaq;
//                
//                //Funcao Heaviside suave
////                if(x>=(x0 - eps) && x<x0) ptil = (eps*gamman*(-(6.+ eps*eps*gamman*gamman)*cos(gamman*x0) - 6.*cos(gamman*(x0-eps))) + 12.*(sin(gamman*(eps - x0)) + sin(gamman*x0)))/temp;
////                
////                if(x>=x0 && x<=x1) ptil = ((cos(gamman*x0)-cos(gamman*x1)))/gamman;
////                
////                if(x>x1 && x<=(x1+eps)) ptil = (eps*gamman*(6.+ eps*eps*gamman*gamman)*cos(gamman*x1) +
////                                                6.*eps*gamman*cos(gamman*(eps+x1)) + 12.*sin(gamman*x1) - 12.*sin(gamman*(eps+x1)))/temp;
////                ptil = (gammaq*ptil*Ft)/gammanq;
//                
//                ptil = (gammaq*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gammanq);
//                util = (gamman*ptil)/gammanq;
//                wtil = (gamman*gamman*ptil)/(gammanq*gammaq);
//                
//                sump += ptil*sin(gamman*x)*sin(gammaq*y);
//                sumux += util*cos(gamman*x)*sin(gammaq*y);
//                sumuy += wtil*sin(gamman*x)*cos(gammaq*y);
//                sumVDx += ptil*gamman*cos(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
//                sumVDy += ptil*gammaq*sin(gamman*x)*cos(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
//            }
//         }
//    }
//    
//    REAL termo1=0., termo2=0.;
//    termo1 = (1.+m);
//    termo2 = 4./(Lx*Ly);
//    
//    if(fabs(y)==0){
//        sump = 2.*sump;
//    }else{
//        sump = termo2*sump;
//    }
//    sumux = -termo1*termo2*sumux;
//    sumuy = termo1*termo2*sumuy;
//    sumVDx = -termo2*sumVDx;
//    sumVDy = -termo2*sumVDy;
//    
//    if(fdimensionless == false){
//        sol[0] = (lamb+2.*mu)*sump;
//        sol[1] = Lref*sumux;
//        sol[2] = Lref*sumuy;
//        flux(0,0) = kovervisc*(lamb+2.*mu)*sumVDx; 
//        flux(1,0) = kovervisc*(lamb+2.*mu)*sumVDy; 
//    }
//    else{
//        sol[0] = sump;
//        sol[1] = sumux;
//        sol[2] = sumuy;
//        flux(0,0) = sumVDx;
//        flux(1,0) = sumVDy;
//    }
//}


void BarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
    
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL eps = 0.0;
    REAL x0 = 0. + eps;
    REAL x1 = Lx - eps;
    REAL m = 1.5;
    REAL tp = ftimeatual;
    REAL Ft = sin(tp);
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
    
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.5*(1.-1./m);
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(5, 0.);
	
    flux(0,0)=0.;
    flux(1,0)=0.;
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil1 = 0., ptil2 = 0.,util = 0.,wtil1 =0., wtil2 =0.;
    REAL sump=0., sump1=0.,sump2=0.,sumux = 0., sumuy = 0., sumVDx = 0., sumVDy=0.;
    
    for(in=1; in<=300; in++)
    {
        gamman = in*PI;
        wtil1 = (2.*(1.+m)*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gamman*gamman);
        sumuy += wtil1*sin(gamman*x);
        
        if(fabs(y)==0)
        {
            ptil1 = ((cos(gamman*x0)-cos(gamman*x1))*Ft)/gamman;
            sump1 += 2.*ptil1*sin(gamman*x);
//            ForcingBCPressao(ptx,sol);
//            sump1 = sol[0];
            sumVDx += 2.*(cos(gamman*x0)-cos(gamman*x1))*Ft*cos(gamman*x);
            sumVDy = 0.;
        }
        
        for(jq=1; jq<=300; jq++)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            
            ptil2 = (4.*gammaq*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gammanq);
            util = ((-1.-m)*gamman*ptil2)/gammanq;
            wtil2 = ((1.+m)*gamman*gamman*ptil2)/(gammanq*gammaq);
            
            if(fabs(y)>0){
                sump2 += ptil2*sin(gamman*x)*sin(gammaq*y);
                sumVDx += ptil2*gamman*cos(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
                sumVDy += ptil2*gammaq*sin(gamman*x)*cos(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
            }
            sumux += util*cos(gamman*x)*sin(gammaq*y);
            sumuy += wtil2*sin(gamman*x)*cos(gammaq*y);
        }
    }
    
    REAL termo1=0.;
    termo1 = 1./(Lx*Ly);
    
    if(fabs(y)==0){
        sump = sump1;
    }else sump = sump2;
    
    sumux = termo1*sumux;
    sumuy = termo1*sumuy;
    sumVDx = -termo1*sumVDx;
    sumVDy = -termo1*sumVDy;
    
    if(fdimensionless == false){
        sol[0] = (lamb+2.*mu)*sump;
        sol[1] = Lref*sumux;
        sol[2] = Lref*sumuy;
        flux(0,0) = kovervisc*(lamb+2.*mu)*sumVDx;
        flux(1,0) = kovervisc*(lamb+2.*mu)*sumVDy; 
    }
    else{
        sol[0] = sump;
        sol[1] = sumux;
        sol[2] = sumuy;
        flux(0,0) = sumVDx;
        flux(1,0) = sumVDy;
    }
}


void SolUBarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv){
    
    
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL eps = 0.0;
    REAL x0 = 0. + eps;
    REAL x1 = Lx - eps;
    REAL m = 1.5;
    REAL tp = ftimeatual;
    REAL Ft = sin(tp);
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
    
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.5*(1.-1./m);
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(2, 0.);
    deriv.Redim(2, 2);
	
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil = 0., util = 0., wtil1 =0., wtil2 =0.;
    REAL sumux = 0., sumuy = 0.;
    REAL sumDux = 0., sumDuy = 0.,sumDuxy = 0.,sumDuyx = 0.;
    
    for(in=1; in<=40; in++)
    {
        gamman = in*PI;
        wtil1 = (2.*(1.+ m)*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gamman*gamman);
        sumuy += wtil1*sin(gamman*x);
        sumDuyx += gamman*wtil1*cos(gamman*x);
            
        for(jq=1; jq<=40; jq++)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            
            ptil = (4.*gammaq*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gammanq);
            util = ((-1.-m)*gamman*ptil)/gammanq;
            wtil2 = ((1.+m)*gamman*gamman*ptil)/(gammanq*gammaq);
            
            sumux += util*cos(gamman*x)*sin(gammaq*y);
            sumuy += wtil2*sin(gamman*x)*cos(gammaq*y);
            
            sumDux += gamman*util*sin(gamman*x)*sin(gammaq*y);//multiplicar por -1
            sumDuxy += gammaq*util*cos(gamman*x)*cos(gammaq*y);
            sumDuy += gammaq*wtil2*sin(gamman*x)*sin(gammaq*y);//multiplicar por -1
            sumDuyx += gamman*wtil2*cos(gamman*x)*cos(gammaq*y);
        }
    }
    
    REAL termo1=0.;
    termo1 = 1./(Lx*Ly);
    
    sumux = termo1*sumux;
    sumuy = termo1*sumuy;
    
    sumDux = -termo1*sumDux;
    sumDuxy = termo1*sumDuxy;
    sumDuy = -termo1*sumDuy;
    sumDuyx = termo1*sumDuyx;
    
    if(fdimensionless == false){
        sol[0] = Lref*sumux;
        sol[1] = Lref*sumuy;
        deriv(0,0) = Lref*sumDux;
        deriv(0,1) = Lref*sumDuxy;
        deriv(1,0) = Lref*sumDuyx;
        deriv(1,1) = Lref*sumDuy;
    }
    else{
        sol[0] = sumux;
        sol[1] = sumuy;
        deriv(0,0) = sumDux;
        deriv(0,1) = sumDuxy;
        deriv(1,0) = sumDuyx;
        deriv(1,1) = sumDuy;
    }
}

void SolPBarryMercerPressureSolution(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &flux){
    
    REAL x = ptx[0];
	REAL y = ptx[1];
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;
    REAL eps = 0.0;
    REAL x0 = 0. + eps;
    REAL x1 = Lx - eps;
    REAL m = 1.5;
    REAL tp = ftimeatual;
    REAL Ft = sin(tp);
    
    REAL perm =0., visc=0., kovervisc=0.;
    
    //dados adimensional
    if(fdimensionless == true)
    {
        perm = 1.;
        visc = 1.;
        kovervisc =perm/visc;
    }
    
    
    //dados com dimensao
    REAL Lref = 0.;
    REAL lamb=0., mu=0.;
    if(fdimensionless == false){
        
        REAL Eyoung = 1.e5;
        REAL poisson = 0.5*(1.-1./m);
        lamb = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
        mu = 0.5*Eyoung/(1.+poisson);
        perm =  1.e-6;
        visc =1.e-4;
        kovervisc =perm/visc;
        Lref = sqrt(Lx*Ly);
        //beta = (lamb+2.*mu)*kovervisc/Lref;
        //tp = tp;
    }
	
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    sol.Resize(5, 0.);
	
    flux(0,0)=0.;
    flux(1,0)=0.;
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
    REAL ptil1=0., ptil2 = 0.;
    REAL sump=0., sump1=0.,sump2=0.,sumVDx = 0., sumVDy=0.;
    
    for(in=1; in<=40; in++)
    {
        gamman = in*PI;
        if(fabs(y)==0)
        {
            ptil1 = ((cos(gamman*x0)-cos(gamman*x1))*Ft)/gamman;
            sump1 += 2.*ptil1*sin(gamman*x);
//            ForcingBCPressao(ptx,sol);
//            sump1 = sol[0];
            sumVDx += 2.*(cos(gamman*x0)-cos(gamman*x1))*Ft*cos(gamman*x);
            sumVDy = 0.;
        }
        
        for(jq=1; jq<=40; jq++)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            
            ptil2 = 4.*(gammaq*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gammanq);
                    
            if(fabs(y)>0){
                sump2 += ptil2*sin(gamman*x)*sin(gammaq*y);
                sumVDx += ptil2*gamman*cos(gamman*x)*sin(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
                sumVDy += ptil2*gammaq*sin(gamman*x)*cos(gammaq*y);//(multiplicado por -1, multipliquei depois da soma)
            }
        }
    }
    
    if(fabs(y)==0){
        sump = sump1;
    }else sump = sump2;
    
    sumVDx = -sumVDx;
    sumVDy = -sumVDy;
    
    if(fdimensionless == false){
        sol[0] = (lamb+2.*mu)*sump;
        //flux(0,0) = kovervisc*(lamb+2.*mu)*sumVDx;
        //flux(1,0) = kovervisc*(lamb+2.*mu)*sumVDy;
    }
    else{
        sol[0] = sump;
        //flux(0,0) = sumVDx;
        //flux(1,0) = sumVDy;
    }
}


REAL HeavisideFunction(REAL val){
    
    if(val>=0.) return 1.;
    else return 0.;
}

void ForcingBCPressao(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    REAL tp = ftimeatual;
    REAL x = pt[0];
    
    REAL BC=0.;
    REAL eps = 0.0;
    REAL x0 = 0. + eps;
    REAL x1 = 1. - eps;
    
    BC = (HeavisideFunction(x-x0) - HeavisideFunction(x-x1));
    
    
//    if(x>=(x0 - eps) && x<=x0) BC = ((eps + 2.*x0 - 2.*x)*(eps - x0 + x)*(eps - x0 + x))/(eps*eps*eps);
//    if(x>=x0 && x<=x1) BC = 1.;
//    if(x>=x1 && x<=(x1+eps)) BC = ((eps + x1 - x)*(eps + x1 - x)*(eps - 2.*x1 + 2.*x))/(eps*eps*eps);
    
    disp[0] = BC*sin(tp);
}

void ForcingBCDeslocamento(const TPZVec<REAL> &pt, TPZVec<STATE> &disp){
    
    REAL x = pt[0];
    
    //dados do problema
	REAL Lx = 1.;

    REAL eps = 0.0;
    REAL x0 = 0. + eps;
    REAL x1 = Lx - eps;
    REAL m = 1.5;
    REAL tp = ftimeatual;
    REAL Ft = sin(tp);
    
    int in, jq;
    REAL PI = atan(1.)*4.;
    
    disp.Resize(4,0.);
	
    
    REAL gamman = 0., gammaq = 0., gammanq = 0.;
	REAL ptil2 = 0., wtil1 =0., wtil2 =0.;
    REAL sumuy = 0., sumvDy = 0., ptil1=0., sump1=0.;
    
    for(in=1; in<=300; in++)
    {
        gamman = in*PI;
        wtil1 = ((cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gamman*gamman);
        sumuy += 2.*wtil1*sin(gamman*x);
        
        ptil1 = ((cos(gamman*x0)-cos(gamman*x1))*Ft)/gamman;
        sump1 += 2.*ptil1*sin(gamman*x);
        
        for(jq=1; jq<=300; jq++)
        {
            gammaq = jq*PI;
            gammanq = gamman*gamman + gammaq*gammaq;
            
            
            ptil2 = (gammaq*(cos(gamman*x0)-cos(gamman*x1))*Ft)/(gamman*gammanq);
            wtil2 = (gamman*gamman*ptil2)/(gammanq*gammaq);
            sumuy += 4.*wtil2*sin(gamman*x);
            sumvDy += ptil2*gammaq*sin(gamman*x)*cos(gammaq);
        }
    }
    
    disp[0] = 0.;//ux
    disp[1] = (1.+m)*sumuy;//uy
    
    REAL BC =(HeavisideFunction(x-x0) - HeavisideFunction(x-x1));
    
//    if(x>=(x0 - eps) && x<=x0) BC = ((eps + 2.*x0 - 2.*x)*(eps - x0 + x)*(eps - x0 + x))/(eps*eps*eps);
//    if(x>=x0 && x<=x1) BC = 1.;
//    if(x>=x1 && x<=(x1+eps)) BC = ((eps + x1 - x)*(eps + x1 - x)*(eps - 2.*x1 + 2.*x))/(eps*eps*eps);
    
    disp[2] = /*sump1;*/BC*sin(tp);//pressao
    
    //disp[3] = -4.*sumvDy;
}

#include "pzgradientreconstruction.h"

//Problema Barry and Mercer: pressure solution
int main_BarryMercerPressureSolution(int argc, char *argv[]){
    
    //    #ifdef LOG4CXX
    //        InitializePZLOG();
    //    #endif
    
    int nthreads = 8;
    
    bool triang = true;
    fdimensionless = true;
    
    ///-------------------
    
    //dados do problema
	REAL Lx = 1.;
    REAL Ly = 1.;

    REAL m = 1.5;//relacao entre os parametros de Lame
    REAL sig0 =0.;
    REAL pini = 0.;
    REAL alpha = 1.0;
    REAL Se = 0.;
    REAL Eyoung = 1.e5;
    REAL poisson = 0.5*(1.-1./m);
    REAL perm = 1.e-6;
    REAL visc = 1.e-4;
   
    REAL rockrho = 0.;
    REAL gravity = 0.;
    REAL fx=0.0;
    REAL fy = gravity*rockrho;
    
    //trabalho com o tempo adimensional
    REAL timeT = /*0.00293215;*/0.00146608/1000;
    REAL lambda = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));//firstlame
    REAL mu = 0.5*Eyoung/(1+poisson);//secondlame
    REAL Se_aux = alpha*alpha/(lambda+2.*mu);
    REAL kovervisc = perm/visc;
    REAL Cf = kovervisc/Se_aux;//fluid diffusivity coeffcient
    REAL beta = Cf/(Ly*Ly);
    timeT = timeT*beta;

    
    //dados adimensional
    if(fdimensionless==true)
    {
        REAL Lref = sqrt(Lx*Ly);
        REAL lambdaD = lambda*Se_aux;
        REAL muD = mu*Se_aux;
        
        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = 0.;
        pini = 0.;
        Lx = Lx/Lref;
        Ly = Ly/Lref;
        perm = 1.;
        visc = 1.;
        fx = 0.;
        fy = 0.;
        Se = 0.;//Incompressible fluid (Se=0)
    }
    
       
    DadosMalhas * mydata = new DadosMalhas();
    mydata->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    
    ofstream saidaerro("Erro.txt");
    for(int p = 1; p<2; p++)
    {
        int pu = p;
        int pq = pu;
        int pp = p-1;
        
        if(triang==false){
            pu = p+1;
            pq = p;
            pp = p;
        }

        int h;
        
        saidaerro<<"\n CALCULO DO ERRO, ELEM. TRIANG., COM ORDEM POLINOMIAL pu = "<< pu << ", pq = "<< pq << " e pp = "<< pp<<endl;
        for (h = 0; h< 5; h++)
        {
            saidaerro<<"\n========= PARA h = "<< h<<"  ============= "<<endl;
            
            // geometric mesh (initial)
            //TPZGeoMesh * gmesh = mydata->GMesh4(Lx, Ly,h,0);//GMesh(triang, Lx, Ly);
            TPZGeoMesh * gmesh = mydata->GMesh(triang, Lx, Ly);
            mydata->UniformRefine(gmesh, h);
            //            std::ofstream malhaGeo("gmesh2D.txt")
            //            gmesh->Print(malhaGeo);
            
            // First computational mesh
            TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu,true,false);
            ofstream arg1("cmesh1.txt");
            cmesh1->Print(arg1);
            
            // second computational mesh
            TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq,true,false);
            ofstream arg2("cmesh2.txt");
            cmesh2->Print(arg2);
            
            // Third computational mesh
            TPZCompMesh * cmesh3=mydata->CMeshPressure(gmesh, pp,triang,true);
            ofstream arg3("cmesh3.txt");
            cmesh3->Print(arg3);
            
            //malha multifisica
            TPZVec<TPZCompMesh *> meshvec(3);
            meshvec[0] = cmesh1;
            meshvec[1] = cmesh2;
            meshvec[2] = cmesh3;
            
            TPZPoroElasticMF2d * mymaterial;
            TPZAutoPointer<TPZFunction<STATE> > BCPress = new TPZDummyFunction<STATE>(ForcingBCPressao, 5);
            TPZAutoPointer<TPZFunction<STATE> > BCTerm = new TPZDummyFunction<STATE>(ForcingBCDeslocamento, 5);
            TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(BarryMercerPressureSolution, 5);
            
            TPZCompMesh * mphysics = mydata->MalhaCompBarryMercerPressureSolution(gmesh,meshvec,mymaterial,BCTerm, solExata);
           // ofstream arg8("mphysic.txt");
            //mphysics->Print(arg8);
            
            //            ofstream arg9("gmesh_final.txt");
            //            gmesh->Print(arg9);
            
            
            int NDeltaT =10000;
            int intervsaidas = NDeltaT/1;
            REAL deltaT=timeT/NDeltaT; //second
            mymaterial->SetTimeStep(deltaT);
            //REAL maxTime = timeT;
            fdeltaT = deltaT;
            
            //mydata->SolveSistBarryMercert(deltaT, maxTime, meshvec, mphysics,1,ftimeatual);
            
            
            ///======== Impimir Erros ========
            TPZAnalysis an(mphysics);
            TPZFMatrix<STATE> Initialsolution = an.Solution();
            
            std::string outputfile;
            outputfile = "TransientSolution";
            
            //Criando matriz de massa (matM)
            TPZAutoPointer <TPZMatrix<STATE> > matM = mydata->MassMatrix(mymaterial, mphysics, nthreads);
            
            //Criando matriz de rigidez (matK) e vetor de carga
            TPZFMatrix<STATE> matK;
            TPZFMatrix<STATE> fvec;
            //mydata->StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec);
            
            int nrows;
            nrows = matM->Rows();
            TPZFMatrix<STATE> TotalRhs(nrows,1,0.0);
            TPZFMatrix<STATE> TotalRhstemp(nrows,1,0.0);
            TPZFMatrix<STATE> Lastsolution = Initialsolution;
            
            REAL TimeValue = 0.0;
            int cent = 1;
            TimeValue = cent*deltaT;
            while (cent<NDeltaT+1)//TimeValue < maxTime || TimeValue == maxTime)
            {
                ftimeatual  = TimeValue;
                // This time solution i for Transient Analytic Solution
                matM->Multiply(Lastsolution,TotalRhstemp);
                mydata->StiffMatrixLoadVec(mymaterial, mphysics, an, matK, fvec, nthreads);
                
                TotalRhs = fvec + TotalRhstemp;
                an.Rhs() = TotalRhs;
                an.Solve();
                Lastsolution = an.Solution();
                
                
                if(cent%intervsaidas==0)
                {
                    saidaerro<<"\n========= PARA O PASSO n = "<< cent <<"  E TEMPO tn = "<< TimeValue <<" =========\n"<<endl;
                    
                    std::stringstream outputfiletemp;
                    outputfiletemp << outputfile << ".vtk";
                    std::string plotfile = outputfiletemp.str();
                    mydata->PosProcessMultphysics(meshvec,mphysics,an,plotfile);
                    
                    TPZVec<REAL> erros;

                    saidaerro<<" Erro da simulacao multifisica do deslocamento (u)" <<endl;
                    TPZAnalysis an12(cmesh1);
                    an12.SetExact(*SolUBarryMercerPressureSolution);
                    bool store_errors = false;
                    an12.PostProcessError(erros, store_errors, saidaerro);

                    saidaerro<<"\nErro da simulacao multifisica da pressao (p)" <<endl;
                    TPZAnalysis an32(cmesh3);
                    an32.SetExact(*SolPBarryMercerPressureSolution);
                    an32.PostProcessError(erros, store_errors, saidaerro);
                }
                
                cent++;
                TimeValue = cent*deltaT;
            }
            
            cmesh1->CleanUp();
            cmesh2->CleanUp();
            cmesh3->CleanUp();
            //mphysics->CleanUp();
            delete cmesh1;
            delete cmesh2;
            delete cmesh3;
            //delete mphysics;
            delete gmesh;
            
            
//            TPZGradientReconstruction *gradreconst = new TPZGradientReconstruction(1);
//            gradreconst->UseWeightCoefficients();
//            //gradreconst->UseSlopeLimiter();
//            cmesh3->LoadReferences();
//            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(meshvec, mphysics);
//            TPZFMatrix<REAL> datagradients;
//            gradreconst->GetDataGradient(cmesh3, datagradients);
//            
//            TPZManVector<std::string,10> scalnames(1), vecnames(0);
//            scalnames[0] = "Solution";
//            std::string plotfile4("saidaSolP.vtk");
//            const int dim = 2;
//            int div = 0;
//            TPZAnalysis ann(cmesh3);
//            ann.DefineGraphMesh(dim,scalnames,vecnames,plotfile4);
//            ann.PostProcess(div,dim);

            
        }
    }
    
   
    
    return 0;
}


void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref){
    
    TPZAutoPointer<TPZRefPattern> refp2D = gRefDBase.FindRefPattern("Qua000022224");
    //TPZAutoPointer<TPZRefPattern> refp1D = gRefDBase.FindRefPattern("Lin000022224");
    
    if(!refp2D) DebugStop();
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < nref; r++)
    {
        int nels = gmesh->NElements();
        for(int iel = 0; iel < nels; iel++)
        {
            gel = gmesh->ElementVec()[iel];
            if(!gel) DebugStop();
            if(gel->Dimension()==2)
            {
                gel->SetRefPattern(refp2D);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
            if(gel->Dimension()==1)
            {
                TPZAutoPointer<TPZRefPattern> refp1D = TPZRefPatternTools::PerfectMatchRefPattern(gel);
                if(refp1D)
                {
                    gel->SetRefPattern(refp1D);
                    TPZVec<TPZGeoEl*> sons;
                    gel->Divide(sons);
                }
                else{
                    DebugStop();//nao conseguiu refinar elementos 1D baseados no refp do pai
                }
            }
            
        }
    }
    
    std::ofstream malhaOut("malhaOut.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
}

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area){
    
    TPZAutoPointer<TPZRefPattern> refpOutroLugar = gRefDBase.FindRefPattern("Qua000022224");
    if(!refpOutroLugar) DebugStop();
    
    int64_t iniEl = 0;
    TPZVec<REAL> qsi(2,0.);
    
    TPZGeoEl * gel = NULL;
    for(int r = 0; r < nref; r++)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl,2);
        if(!gel) DebugStop();
        if(gel->Dimension()==2)
        {
            gel->SetRefPattern(refpOutroLugar);
            TPZVec<TPZGeoEl*> sons;
            gel->Divide(sons);
            for(int edg = gel->NNodes(); edg < gel->NSides(); edg++)
            {
                TPZGeoElSide gelS(gel,edg);
                if(gelS.Dimension() > 1)
                {
                    break;
                }
                TPZGeoElSide neighS(gelS.Neighbour());
                while(neighS != gelS)
                {
                    if(neighS.Element()->Dimension() == 1)
                    {
                        TPZAutoPointer<TPZRefPattern> refp = gel->GetRefPattern();
                        TPZAutoPointer<TPZRefPattern> refp2 = refp->SideRefPattern(gelS.Side());
                        if(refp2)
                        {
                            TPZGeoEl * gel1D = neighS.Element();
                            gel1D->SetRefPattern(refp2);
                            TPZVec<TPZGeoEl*> sons2;
                            gel1D->Divide(sons2);
                        }
                    }
                    neighS = neighS.Neighbour();
                }
            }
        }
    }
    
    if(changeMatId==true)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl,2);
        if(!gel) DebugStop();
        gel->SetMaterialId(newmatId);
        Area = gel->Volume();
    }
    
    std::ofstream malhaOut("malhaOut2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
}

void ErrorHDiv2(TPZCompMesh *hdivmesh, std::ostream &out,void (*fp)(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv))
{
    int64_t nel = hdivmesh->NElements();
    int dim = hdivmesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = hdivmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        cel->EvaluateError(fp, elerror, false);
        int nerr = elerror.size();
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
        
    }
    out << "Errors associated with HDiv space\n";
    out << "L2 Norm for flux = "    << sqrt(globerrors[1]) << endl;
    out << "L2 Norm for divergence = "    << sqrt(globerrors[2])  <<endl;
    out << "Hdiv Norm for flux = "    << sqrt(globerrors[3])  <<endl;
    
}

void ErrorH1(TPZCompMesh *l2mesh, std::ostream &out,void (*fp)(const TPZVec<REAL> &ptx, TPZVec<STATE> &sol, TPZFMatrix<STATE> &deriv))
{
    int64_t nel = l2mesh->NElements();
    int dim = l2mesh->Dimension();
    TPZManVector<REAL,10> globerrors(10,0.);
    for (int64_t el=0; el<nel; el++) {
        TPZCompEl *cel = l2mesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->Dimension() != dim) {
            continue;
        }
        TPZManVector<REAL,10> elerror(10,0.);
        elerror.Fill(0.);
        cel->EvaluateError(fp, elerror, false);
        
        int nerr = elerror.size();
        //globerrors.resize(nerr);
//#ifdef LOG4CXX
//        if (logger->isDebugEnabled()) {
//            std::stringstream sout;
//            sout << "L2 Error sq of element " << el << elerror[0]*elerror[0];
//            LOGPZ_DEBUG(logger, sout.str())
//        }
//#endif
        for (int i=0; i<nerr; i++) {
            globerrors[i] += elerror[i]*elerror[i];
        }
    }
    out << "\n";
    out << "Errors associated with L2 or H1 space\n";
    out << "\nH1 Norm = "    << sqrt(globerrors[0]) << endl;
    out << "\nL2 Norm = "    << sqrt(globerrors[1]) << endl;
    out << "\nSemi H1 Norm = "    << sqrt(globerrors[2]) << endl;
    out << "\n=============================\n"<<endl;
}


