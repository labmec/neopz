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


void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
void TerzaghiProblem1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux);
void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref);
void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area);

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.porolasticmf2d.data"));
#endif

//murad e Loula
int main_murad(int argc, char *argv[])
{
    
#ifdef LOG4CXX
	std::string logs("../logporoelastc2d.cfg");
	InitializePZLOG("../logporoelastc2d.cfg");
#endif
    
    bool triang=true;
    bool dimensionless = true;
    
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
    REAL Ly = 1.;
    REAL Lx = 1.;
    REAL timeT = 10.;
    
    
    //Incompressible fluid (Se=0)
    if(dimensionless==true)
    {
        REAL lambda = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));//firstlame
        REAL mu = 0.5*Eyoung/(1+poisson);//secondlame
        REAL pref = pini;
        REAL Lref = Ly;
        REAL Se_aux = alpha*alpha/(lambda+2.*mu);
        REAL kovervisc = perm/visc;
        REAL Cf = kovervisc/Se_aux;//fluid diffusivity coeffcient
        REAL lambdaD = lambda*Se_aux;
        REAL muD = mu*Se_aux;

        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = sig0/pref;
        pini =pini/pref;
        timeT = timeT*Cf/(Ly*Ly);
        Lx = Lx/Lref;
        Ly = Ly/Lref;
        perm = 1.;
        visc = 1.;
        fx = fx*Lref/pref;
        fy = fy*Lref/pref;
    }
    
    DadosMalhas * mydata = new DadosMalhas();
    mydata->SetParameters(Eyoung, poisson, alpha, Se, perm, visc, fx, fy, sig0);
    
    int pu = 3;
    int pq = 3;
    int pp;
    if(triang==true){
        pp = pq-1;
    }else{
        pq=2;
        pp = pq;
    }
    
	//primeira malha
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = mydata->GMesh(triang,Lx,Ly);
    //mydata->RefiningNearLine(2, gmesh, 3);
    //TPZGeoMesh * gmesh = mydata->GMesh2(LxD,LyD);
    //TPZGeoMesh * gmesh = mydata->GMesh3(LxD,LyD);
    ofstream arg1("gmesh_inicial.txt");
    gmesh->Print(arg1);
    
    // mydata->UniformRefine(gmesh, 2);
    
    
    // First computational mesh
    TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu);
    ofstream arg2("cmesh1_inicial.txt");
    cmesh1->Print(arg2);
    
    // second computational mesh
	TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq);
    ofstream arg3("cmesh2_inicial.txt");
    cmesh2->Print(arg3);
    
    
	// Third computational mesh
	TPZCompMesh * cmesh3 = mydata->CMeshPressure(gmesh, pp,triang);
    ofstream arg4("cmesh3_inicial.txt");
    cmesh3->Print(arg4);
    
    
    // Cleaning reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,4,false);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
    ofstream arg5("cmesh1_final.txt");
    cmesh1->Print(arg5);
	
	
	// Cleaning reference to cmesh2
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,4,false);
   	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
    ofstream arg6("cmesh2_final.txt");
    cmesh2->Print(arg6);
    
    // Cleaning reference to cmesh3
	gmesh->ResetReference();
	cmesh3->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,4,true);
    
	cmesh3->AdjustBoundaryElements();
	cmesh3->CleanUpUnconnectedNodes();
    ofstream arg7("cmesh3_final.txt");
    cmesh3->Print(arg7);
    
    
    //	Set initial conditions for pressure
    TPZAnalysis an3(cmesh3);
	mydata->SolveSist(an3, cmesh3);
	int nrs = an3.Solution().Rows();
    TPZVec<REAL> solini(nrs,pini);
    //	cmesh3->Solution() = solucao1;
    //    cmesh3->LoadSolution(solucao1);
    
    TPZCompMesh  * cmeshL2 = mydata->CMeshPressureL2(gmesh, pp, solini,triang);
    TPZAnalysis anL2(cmeshL2);
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
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(SolucaoExata1D);
    TPZCompMesh * mphysics = mydata->MalhaCompMultphysics(gmesh,meshvec,mymaterial,solExata);
    ofstream arg8("mphysic.txt");
	mphysics->Print(arg8);
    
    REAL deltaT=timeT/100; //secondpra
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = timeT;
    
    
    mydata->SolveSistTransient(deltaT, maxTime, mymaterial, meshvec, mphysics,1,ftimeatual);
    
    return 0;
}

void SolucaoExata1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
    
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
	
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
    flux(0,0)=0.;
    flux(1,0)=0.;
	
	REAL tD = tp;//(lamb+2.*mi)*perm*tp/(visc*H*H);
	REAL xD = fabs(1.-x)/H;
	for (in =0; in<1000; in++) {
		
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
        sol[1] = (-1.)*uD;
        sol[2] = sigD;
        
        flux(1,0) = (1.)*VDy;
    }else{
    
        sol[0] = pD*pini;
        sol[1] = (-1.)*uD*(pini*H)/(lamb+2.*mi);
        sol[2] = (sigD)*pini;
    
        flux(1,0) = (1.)*VDy*pini*(perm/visc);
    }
}


//problema de Terzaghi com Se!=0
int main_Terz(int argc, char *argv[]){
#ifdef LOG4CXX
	std::string logs("../logporoelastc2d.cfg");
	InitializePZLOG("../logporoelastc2d.cfg");
#endif
    
    bool triang=false;
    bool dimensionless=false;
    
    REAL Eyoung = 1.e5;
    REAL poisson = 0.2;
    REAL firstlame = (Eyoung*poisson)/((1.+poisson)*(1.-2.*poisson));
    REAL secondlame = 0.5*Eyoung/(1.+poisson);
    REAL alpha = 1.0;
    REAL Se = 1.e-1;
    
    REAL K = firstlame + (2./3.)*secondlame;
    REAL Ku = K + alpha*alpha/Se;
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
    REAL timeT = 1000;
    
    
    //Incompressible fluid (Se=0)
    if(dimensionless==true)
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
        REAL Cf = kovervisc/Se_aux;//fluid diffusivity coeffcient
        
        //adimemsionalizando
        Eyoung = muD*(3.*lambdaD+2.*muD)/(lambdaD+muD);
        poisson = 0.5*lambdaD/(lambdaD+muD);
        sig0 = sig0/pref;
        pini =pini/pref;
        timeT = timeT*Cf/(Ly*Ly);
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
    
    int pu = 3;
    int pq = 3;
    int pp;
    if(triang==true){
        pp = pq-1;
    }else{
        pq=2;
        pp = pq;
    }
    
	// geometric mesh (initial)
	TPZGeoMesh * gmesh = mydata->GMesh(triang,Lx,Ly);
    //mydata->RefiningNearLine(2, gmesh, 4);

    
    // First computational mesh
    TPZCompMesh * cmesh1 = mydata->MalhaCompElast(gmesh,pu);
    
    // second computational mesh
	TPZCompMesh * cmesh2= mydata->CMeshFlux(gmesh, pq);
    
	// Third computational mesh
    TPZCompMesh * cmesh3=mydata->CMeshPressure(gmesh, pp,triang);
    
    // Cleaning reference of the geometric mesh to cmesh1
	gmesh->ResetReference();
	cmesh1->LoadReferences();
    TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh1,4,false);
	cmesh1->AdjustBoundaryElements();
	cmesh1->CleanUpUnconnectedNodes();
	
	
	// Cleaning reference to cmesh2
	gmesh->ResetReference();
	cmesh2->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh2,4,false);
   	cmesh2->AdjustBoundaryElements();
	cmesh2->CleanUpUnconnectedNodes();
    
    // Cleaning reference to cmesh3
	gmesh->ResetReference();
	cmesh3->LoadReferences();
	TPZBuildMultiphysicsMesh::UniformRefineCompMesh(cmesh3,4,true);
    
	cmesh3->AdjustBoundaryElements();
	cmesh3->CleanUpUnconnectedNodes();
    
    
    //	Set initial conditions for pressure
    TPZAnalysis an3(cmesh3);
	int nrs = an3.Solution().Rows();
    TPZVec<REAL> solini(nrs,pini);
    
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
    TPZAutoPointer<TPZFunction<STATE> > solExata = new TPZDummyFunction<STATE>(TerzaghiProblem1D);
    TPZCompMesh * mphysics = mydata->MalhaCompTerzaghi(gmesh,meshvec,mymaterial,solExata);
    
    REAL deltaT=timeT/1000000; //second
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = timeT;
    
    
    mydata->SolveSistTransient(deltaT, maxTime, mymaterial, meshvec, mphysics,100,ftimeatual);
    
    return 0;
}

void TerzaghiProblem1D(const TPZVec<REAL> &ptx, TPZVec<REAL> &sol, TPZFMatrix<REAL> &flux){
    
    bool dimensionless=false;
    
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
	
	sol[0]=0.;
	sol[1]=0.;
	sol[2]=0.;
    flux(0,0)=0.;
    flux(1,0)=0.;
	
    REAL aux1= (alpha*F)/(Se*(Ku + (4./3.)*mu));
    REAL aux2=0.;
    REAL aux3=0.;
    
    if(dimensionless==true)
    {
        for (in =0; in<1000; in++)
        {
            aux2 = 4./(PI*(2.*in + 1.));
            aux3 = Pi*(2.*in + 1.)/2.;
            
            sump += aux2*sin(aux3*y)*exp(-aux3*aux3*tp);
            
            sumuyH += 0.5*aux2*aux2*exp(-aux3*aux3*tp);
            
            sumVDy += cos(aux3*y)*exp(-aux3*aux3*tp);
        }
        
        p = sump;
        uyH = SInfty + (S0 - SInfty)*sumuyH/aux1;
        VDy = (-2.)*sumVDy;
        
        sol[0] = p;
        sol[1] = -uyH;
        flux(1,0) = VDy;
        
    }else{
        for (in =0; in<1000; in++)
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
        sol[1] = -uyH;
        flux(1,0) = -VDy;
    }
}


int main(int argc, char *argv[]){

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
        gRefDBase.InsertRefPattern(refp);
        if(!refp)
        {
            DebugStop();
        }
    }
    
    
    REAL Lx = 1.;
    REAL Ly = 1.;
    DadosMalhas * mydata = new DadosMalhas();
    TPZGeoMesh * gmesh = mydata->GMesh(false,Lx,Ly);
    
    RefinamentoPadrao3x3(gmesh,0);
    
    TPZVec<REAL> pt(3);
    pt[0] = Lx/4.;
    pt[1] = Ly/4.;
    pt[2] = 0.;
    int newmatId = mydata->GetIdSourceTerm();//mat id of the source term
    REAL Area;
    RefinamentoPadrao3x3(gmesh, 2,pt, true, newmatId, Area);
    
    std::ofstream malhaGeo("gmesh2D.txt");
    gmesh->Print(malhaGeo);

    
    return 0;
}

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref){
    
   TPZAutoPointer<TPZRefPattern> refpOutroLugar = gRefDBase.FindRefPattern("Qua000022224");
    if(!refpOutroLugar) DebugStop();
  
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
                gel->SetRefPattern(refpOutroLugar);
                TPZVec<TPZGeoEl*> sons;
                gel->Divide(sons);
            }
        }
    }
       
    std::ofstream malhaOut("malhaOut.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);
}

void RefinamentoPadrao3x3(TPZGeoMesh *gmesh, int nref,TPZVec<REAL> pt, bool changeMatId, int newmatId, REAL &Area){

    TPZAutoPointer<TPZRefPattern> refpOutroLugar = gRefDBase.FindRefPattern("Qua000022224");
    if(!refpOutroLugar) DebugStop();
       
    int iniEl = 0;
    TPZVec<REAL> qsi(2,0.);

    TPZGeoEl * gel = NULL;
    for(int r = 0; r < nref; r++)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl);
        if(!gel) DebugStop();
        if(gel->Dimension()==2)
        {
            gel->SetRefPattern(refpOutroLugar);
            TPZVec<TPZGeoEl*> sons;
            gel->Divide(sons);
        }
    }
    
    if(changeMatId==true)
    {
        gel = gmesh->FindElement(pt, qsi, iniEl);
        if(!gel) DebugStop();
        gel->SetMaterialId(newmatId);
        Area = gel->Volume();
    }
                
    std::ofstream malhaOut("malhaOut2.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh, malhaOut, true);

    
}
