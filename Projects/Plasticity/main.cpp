
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

//#include "TPZPlasticityTest.h"
#include <iostream>
#include <cstdlib>
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "TPZTensor.h"
#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZYCVonMises.h"
#include "TPZVonMises.h"
#include "pzfstrmatrix.h"
#include "pzbndmat.h"
#include "pzgeoquad.h"
#include "TPZGeoCube.h"
#include "pzgeotetrahedra.h"
#include "pzgeopyramid.h"
#include "tpzgeoelrefpattern.h"
#include "pzbndcond.h"
#include "pzstepsolver.h"
#include "TPZTensor.h"
#include "TPZYCMohrCoulomb.h"
#include "TPZMohrCoulomb.h"
#include "TPZDruckerPrager.h"
#include "GeoMeshClass.h"
#include "pzelastoplastic2D.h"
#include <pzmathyperelastic.h>
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"


using namespace pzshape; // needed for TPZShapeCube and related classes

static void calcVonMisesBar();
void MohrCoulombTestX();
void cmesh(TPZCompMesh *CMesh, TPZMaterial * mat);//,REAL theta,int axes);
void taludecmesh(TPZCompMesh *CMesh, TPZMaterial * mat);
void calctalude();
TPZGeoMesh * barmesh(int h);
void SolverSetUp(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SetUPPostProcessVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);
void ManageIterativeProcessPesoProprio(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);





#include "pzlog.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.main"));
#endif


int main2()
{
  
    
//    TPZSandlerDimaggio sandler;
//    sandler.McCormicRanchSand(sandler);
//    ofstream outfiletxty("MohrNetoy.txt");
//        ofstream outfiletxtS("SandlerApplyLoad.txt");
//    TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
//    
//    
//    sandler.fIntegrTol = 1.e-12;
//    sandler.SetIntegrTol(1.e-4);
//    
//    deltaeps.XX()= -0.0013;
//    deltaeps.YY()=0;
//    deltaeps.ZZ()=0;
//    eps=deltaeps;
//    
//    
//    for(int i=0;i<100;i++)
//    {
//        
//        sandler.ApplyStrainComputeSigma(eps, sigma);
//       // sandler.ApplyLoad(sigma, eps);
//        
//        TPZPlasticState<REAL> state = sandler.GetState();
//        //state.Print(cout);
//        if(i==49)
//        {
//            deltaeps*=-1;
//        }
//        
//        outfiletxty << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
//        eps+=deltaeps;
//    }
//    
//    TPZSandlerDimaggio sandler2;
//    sandler2.McCormicRanchSand(sandler2);
//    sandler2.SetIntegrTol(0.0000001);
//    deltaeps.XX()= 1.5*-0.065;
//    deltaeps.YY()=0;
//    deltaeps.ZZ()=0;
//    eps=deltaeps;
//    sandler2.ApplyStrainComputeSigma(eps, sigma);
//    //outfiletxty << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
//
//    ///TESTE DE COMPRESSAO SIMPLES
//    deltasigma.XX()=-0.8;
//    deltasigma.YY()=deltasigma.XX()*0.6;
//    deltasigma.ZZ()=deltasigma.YY();
//    sigma=deltasigma;
//    
//    TPZSandlerDimaggio sandler3;
//    sandler3.McCormicRanchSand(sandler3);
//    
//    for(int i=0;i<1;i++)
//    {
//        sandler3.ApplyLoad(sigma,eps);
//        outfiletxtS << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
//         sigma+=deltasigma;
//    }
    
    /*
    
    InitializePZLOG();
      ofstream outfiletxt("MohrNetox.txt");
    ofstream outfiletxty("MohrNetoy.txt");
    //   MohrCoulombTestX();
    //   calctalude();
    //   calcVonMisesBar();
    TPZMohrCoulombNeto toto;
    TPZTensor<REAL> epstotal,sigma,epsplastic;
    TPZFNMatrix<36,REAL> tangent(6,6);
    REAL deltaeps =0.00005;
    epstotal.XX() = deltaeps;
    TPZFMatrix<REAL> curve(7,2);

    for(int i =1;i<360;i++)
    {
        cout << "\n i = "<<i<<endl;
        TPZMohrCoulombNeto::TComputeSequence memory;
        memory = toto.ComputeSigma<REAL>(epstotal, sigma);
        toto.CommitDeformation(epstotal,memory);
        epstotal.XX()+=deltaeps;
        //epsplastic=toto.fState.fEpsPlastic;
        outfiletxty << epstotal.XX() << " " << sigma.XX() << "\n";
        if(i==50 || i==80)
        {
            deltaeps*=-1;
        }
        
    }
    
    
    do {
        TPZMohrCoulombNeto::TComputeSequence memory;
        memory = toto.ComputeSigma(epstotal, sigma );
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Deformation tensor ";
            epstotal.Print(sout);
            sout << "Stress tensor " << sigma;
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        
        toto.ComputeSigmaTangent(epstotal, sigma, tangent, memory);
        toto.CommitDeformation(epstotal, memory);
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "Total deformation tensor ";
            epstotal.Print(sout);
            sout << "Stress tensor " << sigma;
            tangent.Print("Tangent stress",sout);
            toto.Print(sout);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        epstotal.XX() += deltaeps;
      
        outfiletxt << fabs(epstotal.XX()) << " " << fabs(sigma.XX()) << "\n";
        
        
    
    } while (epstotal.XX() < 0.0035);

    */
    return 0;
}


static void calcVonMisesBar()
{
    TPZGeoMesh *barmesh1;
    barmesh1 = barmesh(0);
    ofstream arg("barmesh.txt");
    barmesh1->Print(arg);
    
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(barmesh1);
	
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
    TPZVonMises VM;
    TPZVonMises::Steel(VM);
    
	TPZMatElastoPlastic<TPZVonMises> PlasticVonM(1);
	PlasticVonM.SetPlasticity(VM);
    
    TPZMaterial *plastic(&PlasticVonM);
    compmesh1->InsertMaterialObject(plastic);
    
//    REAL theta = 45;
//    int axes=2;
    //cmesh(compmesh1,plastic,theta,axes);
    // FIX ME
    DebugStop();
    
    ofstream arg2("barcmesh.txt");
    compmesh1->Print(arg2);
    
	
	TPZElastoPlasticAnalysis EPAnalysis(compmesh1,cout);
    SolverSetUp(EPAnalysis,compmesh1);
	
	TPZPostProcAnalysis PPAnalysis(compmesh1);
	TPZFStructMatrix structmatrix(PPAnalysis. Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
    
    PPAnalysis.SetStructuralMatrix(structmatrix);
    TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
    
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =3;
    
    
	
    //	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"barraerick4.vtk");
    //  PPAnalysis.PostProcess(0/*pOrder*/);
	
    TPZFMatrix<REAL> BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix<REAL> val1(3,1,0.);TPZFMatrix<REAL> val2(3,1,0.);TPZFMatrix<REAL> BeginForce(3,1,0.);TPZFMatrix<REAL> EndForce(3,1,0.);
	
	int nsteps,taxa,nnewton;
	nnewton = 30;
	nsteps =1;
	taxa = 1;
	REAL beginforce = 209.9999;
	REAL endforce = 209.9999;
    int bc =-2;
    
    ManageIterativeProcess(EPAnalysis ,beginforce,endforce,bc, nsteps,&PPAnalysis);
    
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"barraerick8.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();

}


void calctalude()
{
    TPZGeoMesh *TaludeGmesh;
    TaludeGmesh = GeoMeshClass::Talude();
    ofstream arg("TaludeGmesh.txt");
    TaludeGmesh->Print(arg);
    
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(TaludeGmesh);
	
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
    TPZMohrCoulomb MC;
    TPZMohrCoulomb::TaludeMaterial(MC);
    
    TPZDruckerPrager DP;
    TPZDruckerPrager::TaludeMaterial(DP,0);
    
    TPZMatElastoPlastic2D<TPZDruckerPrager> PlasticDP(1,1);
    PlasticDP.SetPlasticity(DP);
    PlasticDP.SetBulkDensity(1.);//kN/mˆ3
    
    
	TPZMatElastoPlastic2D<TPZMohrCoulomb> PlasticMohrC(1,1);
	PlasticMohrC.SetPlasticity(MC);
    PlasticMohrC.SetBulkDensity(1.);//kN/mˆ3
    
   // TPZMaterial *plastic(&PlasticMohrC);
    TPZMaterial *plasticDP(&PlasticDP);
    //compmesh1->InsertMaterialObject(plastic);
    compmesh1->InsertMaterialObject(plasticDP);
    
    taludecmesh(compmesh1,plasticDP);
    
    ofstream arg2("barcmesh.txt");
    compmesh1->Print(arg2);
    
	
	TPZElastoPlasticAnalysis EPAnalysis(compmesh1,cout);
    SolverSetUp(EPAnalysis,compmesh1);
	
	TPZPostProcAnalysis PPAnalysis(compmesh1);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
    
    PPAnalysis.SetStructuralMatrix(structmatrix);
    TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables(PostProcVars,scalNames,vecNames);
    
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
    
    
	
    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"talude.vtk");
    PPAnalysis.PostProcess(0/*pOrder*/);
	
    TPZFMatrix<REAL> BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix<REAL> val1(3,1,0.);TPZFMatrix<REAL> val2(3,1,0.);TPZFMatrix<REAL> BeginForce(3,1,0.);TPZFMatrix<REAL> EndForce(3,1,0.);
	
	int nsteps,taxa,nnewton;
	nnewton = 5;
	nsteps =5;
	taxa = 1;
	REAL beginforce = 0.1;
	REAL endforce = 0.2;
    int bc =1;
    
    ManageIterativeProcessPesoProprio(EPAnalysis ,beginforce,endforce,bc, nsteps,&PPAnalysis);
    
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"talude.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	

}



void taludecmesh(TPZCompMesh *CMesh, TPZMaterial * mat)
{
    
    TPZFMatrix<REAL> k1(3,3,0.);
    TPZFMatrix<REAL> f1(3,1,0.);
    TPZMaterial *bc1 = mat->CreateBC(mat,-1,0,k1,f1);
    CMesh->InsertMaterialObject(bc1);
    
    TPZFMatrix<REAL> k2(3,3,0.);
    TPZFMatrix<REAL> f2(3,1,0.);
    f2(0,0)=1.;
    TPZMaterial *bc2 = mat->CreateBC(mat,-2,3,k2,f2);
    CMesh->InsertMaterialObject(bc2);
    
    TPZFMatrix<REAL> k3(3,3,0.);
    TPZFMatrix<REAL> f3(3,1,0.);
    f3(0,0)=1.;
    TPZMaterial *bc3 = mat->CreateBC(mat,-3,3,k3,f3);
    CMesh->InsertMaterialObject(bc3);
    
    
    CMesh->AutoBuild();
    
    
}


void cmesh(TPZCompMesh *CMesh, TPZMaterial * mat)//,REAL theta,int axes)
{
    
    TPZFMatrix<REAL> k1(3,3,0.);
    TPZFMatrix<REAL> f1(3,1,0.);
    TPZMaterial *bc1 = mat->CreateBC(mat,-1,0,k1,f1);
    CMesh->InsertMaterialObject(bc1);
    TPZFMatrix<REAL> k2(3,3,0.);
    TPZFMatrix<REAL> f2(3,1,0.);
    f2(0,0)=1.;

    TPZMaterial * bc2 = mat->CreateBC(mat,-2,1,k2,f2);
    CMesh->InsertMaterialObject(bc2);
    
    CMesh->AutoBuild();
    
}


TPZGeoMesh * barmesh(int h)
{
    TPZGeoMesh * gMesh = new TPZGeoMesh;
    TPZVec<TPZGeoNode> nodes(8);
    gMesh->NodeVec().Resize(8);
    
    nodes[0].SetNodeId(0);
    nodes[0].SetCoord(0,0.);
    nodes[0].SetCoord(1,0.);
    nodes[0].SetCoord(2,0.);
    gMesh->NodeVec()[0] = nodes[0];
    
    nodes[1].SetNodeId(1);
    nodes[1].SetCoord(0,0.);
    nodes[1].SetCoord(1,0.);
    nodes[1].SetCoord(2,-100.);
    gMesh->NodeVec()[1] = nodes[1];
    
    nodes[2].SetNodeId(2);
    nodes[2].SetCoord(0,0.);
    nodes[2].SetCoord(1,100.);
    nodes[2].SetCoord(2,-100.);
    gMesh->NodeVec()[2] = nodes[2];
    
    nodes[3].SetNodeId(3);
    nodes[3].SetCoord(0,0.);
    nodes[3].SetCoord(1,100.);
    nodes[3].SetCoord(2,0.);
    gMesh->NodeVec()[3] = nodes[3];
    
    nodes[4].SetNodeId(4);
    nodes[4].SetCoord(0,1000.);
    nodes[4].SetCoord(1,0.);
    nodes[4].SetCoord(2,0.);
    gMesh->NodeVec()[4] = nodes[4];
    
    nodes[5].SetNodeId(5);
    nodes[5].SetCoord(0,1000.);
    nodes[5].SetCoord(1,0.);
    nodes[5].SetCoord(2,-100.);
    gMesh->NodeVec()[5] = nodes[5];
    
    nodes[6].SetNodeId(6);
    nodes[6].SetCoord(0,1000.);
    nodes[6].SetCoord(1,100.);
    nodes[6].SetCoord(2,-100.);
    gMesh->NodeVec()[6] = nodes[6];
    
    nodes[7].SetNodeId(7);
    nodes[7].SetCoord(0,1000.);
    nodes[7].SetCoord(1,100.);
    nodes[7].SetCoord(2,0.);
    gMesh->NodeVec()[7] = nodes[7];
    
    
    
    TPZVec<long> barra(8),quad1(4),quad2(4);
    
    barra[0]=0;
    barra[1]=1;
    barra[2]=2;
    barra[3]=3;
    barra[4]=4;
    barra[5]=5;
    barra[6]=6;
    barra[7]=7;
    
    quad1[0]=0;
    quad1[1]=1;
    quad1[2]=2;
    quad1[3]=3;
    
    quad2[0]=4;
    quad2[1]=5;
    quad2[2]=6;
    quad2[3]=7;
    
    new TPZGeoElRefPattern< pzgeom::TPZGeoCube > (1,barra,1,*gMesh);
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (2,quad1,-1,*gMesh);
    new TPZGeoElRefPattern< pzgeom::TPZGeoQuad > (3,quad2,-2,*gMesh);
    
    
    gMesh->BuildConnectivity();
    
    //    for(int ref = 0; ref <0; ref++)
    //	{
    //		TPZVec<TPZGeoEl *> tatara;
    //		int n = gMesh->NElements();
    //		for(int i = 0; i < n; i++)
    //		{
    //			TPZGeoEl * gel = gMesh->ElementVec()[i];
    //			gel->Divide(tatara);
    //		}
    //	}
    
    //ofstream barmesh("barmesh.vtk");
    return gMesh;
    
}

void SolverSetUp(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    
    //TPZFStructMatrix full(fCmesh)
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
    //  step.SetDirect(ELDLt);
    //  step.SetJacobi(5000, 1.e-12,0);
    step.SetDirect(ECholesky);
    // step.SetDirect(ELU);
	an.SetSolver(step);
	
}

void SetUPPostProcessVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	scalnames.Resize(8);
	scalnames[0] = "Alpha";
	scalnames[1] = "PlasticSteps";
	scalnames[2] = "VolElasticStrain";
	scalnames[3] = "VolPlasticStrain";
	scalnames[4] = "VolTotalStrain";
	scalnames[5] = "I1Stress";
	scalnames[6] = "J2Stress";
	scalnames[7] = "YieldSurface";
    //	scalnames[8] = "EMisesStress";
    
	vecnames.Resize(5);
	vecnames[0] = "Displacement";
	vecnames[1] = "NormalStress";
	vecnames[2] = "ShearStress";
	vecnames[3] = "NormalStrain";
	vecnames[4] = "ShearStrain";
    
 	postprocvars.Resize(scalnames.NElements()+vecnames.NElements());
	int i, k=0;
	for(i = 0; i < scalnames.NElements(); i++)
	{
		postprocvars[k] = scalnames[i];
		k++;
	}
	for(i = 0; i < vecnames.NElements(); i++)
	{
		postprocvars[k] = vecnames[i];
		k++;
	}
}

void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost)
{

    TPZGeoMesh *barmesh1;
    barmesh1 = barmesh(0);
    ofstream arg("barmesh.txt");
    barmesh1->Print(arg);
    
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(barmesh1);
	
    
	//TPZAnalysis::SetAllCreateFunctionsContinuous(compmesh1);
  
    
    int nummat=1;
    REAL e = 5088.; // [MPa]
    REAL mu = 2290.;
    REAL nu = 0.11;
    REAL lambda = 646;
    
    REAL coef1;// =1.202;
    coef1=mu*((1-3*nu)/(1-2*nu));
    REAL coef2;// = -0.057;
    coef2 = mu*((1-nu)/(1-2*nu));
    REAL coef3;// = 0.004;
    coef3 = mu*((2*nu)/(1-2*nu));
    
    TPZMatHyperElastic * mathyper = new TPZMatHyperElastic(nummat,e,mu,nu,lambda, coef1,coef2,coef3);
	TPZMaterial *hyperm(mathyper);
    compmesh1->InsertMaterialObject(hyperm);
    cmesh(compmesh1, hyperm);
    ofstream arg2("hyper.txt");
    compmesh1->Print(arg2);
    
    TPZNonLinearAnalysis an(compmesh1,cout);
    
    TPZSkylineStructMatrix full(compmesh1);
	an.SetStructuralMatrix(full);
    an.Solution().Zero();
    TPZStepSolver<STATE> sol;
    sol.SetDirect(ELDLt);
    an.SetSolver(sol);
    an.IterativeProcess(cout,1e-5,30);
    TPZStack<std::string> scalnames,vecnames;
    vecnames.Push("Solution");
    vecnames.Push("VonMises");
    an.DefineGraphMesh(3, scalnames, vecnames, "hyper.vtk");
    int postprocessresolution = 0;
    an.PostProcess(postprocessresolution);
    
    
	TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    
	if(!pBC)return;
	
    //    analysis.Mesh()->MaterialVec()[1]->SetBulkDensity(1.);
    //mat->
    
    REAL increment =  (valEnd-valBeg)/steps;
    cout<< "\n increment = " <<increment<<endl;
    TPZFMatrix<REAL> mattemp(3,1,0.);
    mattemp(0,0)=valBeg;
    int i;
    
    
    for(i = 0; i <= steps; i++)
	{
        
		
        mattemp(0,0)+=increment;
        cout << "\n mattemp "<<mattemp<<endl;
        
		pBC->Val2() = mattemp;
        cout<< "\n pBC->Val2() = " <<pBC->Val2()<<endl;
        bool linesearch = false;
        bool checkconv = false;
		analysis.IterativeProcess(cout, 1.e-5, 30, linesearch, checkconv);
        
        analysis.AcceptSolution();
        cout << "\nPostSolution "<< endl<< pPost->Solution();
        analysis.TransferSolution(*pPost);
        pPost->PostProcess(0);
        cout << pPost->Solution();
        
	}
}

void ManageIterativeProcessPesoProprio(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost)
{
    
    
	TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
    //TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    TPZMatElastoPlastic<TPZMohrCoulomb> * pMC = dynamic_cast<TPZMatElastoPlastic<TPZMohrCoulomb> *>(mat);
    // TPZMatElastoPlastic<TPZDruckerPrager> * pMC = dynamic_cast<TPZMatElastoPlastic<TPZDruckerPrager> *>(mat);
    
    REAL increment =  (valEnd-valBeg)/steps;
    cout<< "\n increment = " <<increment<<endl;
    TPZFMatrix<REAL> mattemp(3,1,0.);
    mattemp(0,0)=valBeg;
    int i;
    
    
    for(i = 0; i <= steps; i++)
	{
        
		
        mattemp(0,0)+=increment;
        REAL tempmattemp= mattemp(0,0);
        cout << "\n mattemp "<<mattemp<<endl;
        
        //  mat->SetBulkDensity(increment);
        
		pMC->SetBulkDensity(tempmattemp);
        bool linesearch = false;
        bool checkconv = false;
		analysis.IterativeProcess(cout, 1.e-5, 5, linesearch, checkconv);
        
        analysis.AcceptSolution();
        cout << "\nPostSolution "<< endl<< pPost->Solution();
        analysis.TransferSolution(*pPost);
        pPost->PostProcess(0);
        // cout << pPost->Solution();
        
	}
}

#include "TPZYCModifiedMohrCoulomb.h"

void MohrCoulombTestX()
{
	
	
    
    
	ofstream outfiletxt("DPInTriaxialCompression.txt");
	
	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
	
    deltastress.XX() =-1.;
    deltastress.XY() = 0.;
    deltastress.XZ() = 0.;
    deltastress.YY() = -1.;
    deltastress.YZ() = 0.;
    deltastress.ZZ() = -0.1;
    stress = deltastress;
    
    //3 2 1 cai no inner
    //1 2 3 cai no inner
    //2 3 1 cai no inner
    //2 1 6 cai no outer
    //2 1 0 cai MC
	
    //	TPZFNMatrix<6*6> Dep(6,6,0.);
    //    deltastrain.XX() = -0.0002;
    //	deltastrain.XY() = 0.;
    //	deltastrain.XZ() = 0.;
    //	deltastrain.YY() = -0.00000001;
    //	deltastrain.YZ() = 0.;
    //	deltastrain.ZZ() = -0.00000001;
    //	strain=deltastrain;
	
    TPZPlasticStep<TPZYCModifiedMohrCoulomb,TPZThermoForceA,TPZElasticResponse> modifiedmohr;
    
    REAL pi = M_PI;
    REAL cohesion = 11.2033; //yield- coesao inicial correspondeno a fck igual 32 Mpa
    REAL phi =  20./180. * pi; //phi=20
    REAL hardening = 1000.; //Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao
    REAL young = 20000.;
    REAL poisson = 0.;
    modifiedmohr.fYC.SetUp(phi);
    modifiedmohr.fTFA.SetUp(cohesion, hardening);
    modifiedmohr.fER.SetUp(young, poisson);
    
	TPZMohrCoulomb Pstep;
    Pstep.ConventionalConcrete(Pstep);
    
    TPZDruckerPrager DP;
    DP.ConventionalConcrete(DP,0);
    
    int length =430;
	for(int step=0;step<length;step++)
	{
		cout << "\nstep "<< step;
        
		DP.ApplyLoad(stress,strain);
        stress += deltastress;
        //		cout<<  "\nstress " << stress << endl;
        //		cout<<  "\nstrain " << strain << endl;
        //		Pstep.Print(cout);
        cout << "\n strain " << strain << endl;
		TPZVec<REAL> phis(1);
		DP.Phi(strain, phis);
        
		cout << "\nphis " << phis << endl;
 		//cout<<  "\nEigen " << eigenval << endl;
        //if(step==50 || step == 100 ||step == 200 ||step == 300)deltastress*=-1.;
		outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n"; 
		
	}
	
	
}

