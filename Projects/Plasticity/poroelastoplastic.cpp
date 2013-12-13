#include "poroelastoplastic.h"
#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
#include "TPZLadeKim.h"
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
#include "TPZProjectEllipse.h"

#include "WellBoreAnalysis.h"

void VisualizeSandlerDimaggio(std::stringstream &FileName, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD);

void SolverSetUp2(TPZAnalysis &an, TPZCompMesh *fCmesh);
void SetUPPostProcessVariables2(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void ManageIterativeProcess2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);
void ManageIterativeProcessPesoProprio2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost);
void RotationMatrix(TPZFMatrix<REAL> &R, double thetaRad, int axis);
void RotateMatrix(TPZFMatrix<REAL> &Mat, double thetaRad,int rotateaboutaxes);
void RotateMesh(TPZGeoMesh &geomesh, REAL angleDegree,int rotateabout);
void calcSDBar();
TPZGeoMesh * BarMesh(int h);
void Cmesh(TPZCompMesh *CMesh, TPZMaterial * mat,REAL theta,int axes);
void wellcmesh();
void wellboreanalysis();

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("plasticity.poroelastoplastic"));
#endif

using namespace pzshape; // needed for TPZShapeCube and related classes

#include <math.h>

void wellboreanalyis()
{

    TPZGeoMesh *gmesh = new TPZGeoMesh();
    GeoMeshClass::WellBore2d(gmesh);
    ofstream arg("wellgeomeshlog.txt");
    gmesh->Print(arg);
    
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(gmesh);
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> SD;
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::UncDeepSandTest(SD);
    SD.SetResidualTolerance(1.e-10);
    SD.fIntegrTol = 1.;
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
	TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> > PlasticSD(1,1);
    
    TPZTensor<STATE> Confinement;
    Confinement.XX() = -44.3;
    Confinement.YY() = -58.2;
    Confinement.ZZ() = -53.8;
    
    TPZTensor<REAL> initstress,finalstress;
    REAL hydro = Confinement.I1();
    hydro -= SD.fYC.fA*SD.fYC.fR;
    hydro /= 3.;
    finalstress.XX() = hydro;
    finalstress.YY() = hydro;
    finalstress.ZZ() = hydro;
    
    PrepareInitialMat(SD, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = Confinement;
    PrepareInitialMat(SD, initstress, finalstress, 10);
    
    TPZMaterial *plastic(&PlasticSD);
    
    PlasticSD.SetPlasticity(SD);
    compmesh1->InsertMaterialObject(plastic);
    
    CmeshWell(compmesh1,plastic,Confinement,hydro);
  
    

    
    TPZElastoPlasticAnalysis analysis(compmesh1,cout);

    // Eh preciso colocar depois da criacao do objeto analysis porque o sequence number muda na construcao do analyse
    TPZCompMesh locmesh(*compmesh1);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        compmesh1->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    

    SolverSetUp2(analysis,compmesh1);
    
    
    
 
    
    int BCId=-2;
    TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    int steps = 5;
    
    TPZFNMatrix<9,STATE> mattemp(3,3,0.), matinit(3,3,0.),matfinal(3,3,0.), matincrement(3,3,0.);
    matinit = pBC->Val1();
    matfinal(0,0) = -29.3;
    matfinal(1,1) = -29.3;

    matincrement = matfinal-matinit;
    matincrement *= (1./steps);
    mattemp = matinit;
    
    std::string vtkFile = "pocoplastico.vtk";
    TPZPostProcAnalysis ppanalysis(analysis.Mesh());
    TPZFStructMatrix structmatrix(ppanalysis.Mesh());
    structmatrix.SetNumThreads(8);
    ppanalysis.SetStructuralMatrix(structmatrix);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZStack<std::string> PostProcVars, scalNames, vecNames;
    SetUPPostProcessVariables2(PostProcVars,scalNames, vecNames);
    //    
    ppanalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
    //
    ppanalysis.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
    //	
    
    matincrement.Print("Incremento de tensao",std::cout);
    
    int neq = analysis.Mesh()->Solution().Rows();
    TPZFMatrix<STATE> allsol(neq,steps+1,0.);
    
    
    for(int i=0;i<=steps;i++)
    {
        pBC->Val1()=mattemp;
        bool linesearch = false;
        bool checkconv = false;
        analysis.IterativeProcess(cout, 1.e-8, 30, linesearch, checkconv);
        
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
        
        
        ////Post Processing
        TPZVec<std::string> vecnames,scalnames;
        scalnames.Resize(0);
        vecnames.Resize(1);
        vecnames[0]="Displacement";
        // vecnames[1]="NormalStrain";

        TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            allsol(ieq,i) = sol(ieq,0);
        }
    
        analysis.AcceptSolution();
        
        
        analysis.TransferSolution(ppanalysis);
        ppanalysis.PostProcess(0);// pOrder
        

        mattemp += matincrement;
        
    }
    
    locmesh.Solution() = allsol;
    

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (locmesh.MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    TPZFMatrix<STATE> locrhs(neq,steps,0.);
    TPZSkylineStructMatrix skylstr(&locmesh);
    TPZMatrix<STATE> *matrix = skylstr.TPZStructMatrix::CreateAssemble(locrhs, 0);
    delete matrix;
    locmesh.Solution() = analysis.CumulativeSolution();
    pMatWithMem2->SetUpdateMem(false);

/*
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem1 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (analysis.Mesh()->MaterialVec()[1]);
    int elnum = 89;
    TPZCompEl *cel1 = analysis.Mesh()->ElementVec()[elnum];
    TPZCompEl *cel2 = locmesh.ElementVec()[elnum];
    TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    TPZMaterialData data1,data2;
    data1.intPtIndex = 0;
    data2.intPtIndex = 0;
    TPZManVector<REAL,3> qsi(2,0.);
    REAL weight;
    intel1->GetIntegrationRule().Point(0, qsi, weight);
    intel1->InitMaterialData(data1);
    intel2->InitMaterialData(data2);
    intel1->ComputeRequiredData(data1, qsi);
    intel2->ComputeRequiredData(data2, qsi);
    TPZManVector<STATE,4> solout2(1),solout1(1);
    int var = pMatWithMem2->VariableIndex("PlasticSqJ2");
    pMatWithMem2->Solution(data2, var, solout2);
    pMatWithMem1->Solution(data1, var, solout1);
*/    
    {
        std::string vtkFile = "pocoplastico_copy.vtk";
        
        TPZPostProcAnalysis ppanalysis(&locmesh);
        TPZFStructMatrix structmatrix(ppanalysis.Mesh());
        structmatrix.SetNumThreads(8);
        ppanalysis.SetStructuralMatrix(structmatrix);
        
        TPZVec<int> PostProcMatIds(1,1);
        TPZStack<std::string> PostProcVars, scalNames, vecNames;
        SetUPPostProcessVariables2(PostProcVars,scalNames, vecNames);
        //    
        ppanalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
        //	
        analysis.TransferSolution(ppanalysis);
        //
        ppanalysis.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
        //		
        ppanalysis.PostProcess(0/*pOrder*/);
    }


}

void wellelastic()
{
    REAL E=29269.,nu=0.203;
    int planestrain=1;
    //int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress
    TPZElasticityMaterial *elastic = new TPZElasticityMaterial(1,E,nu,0,0,planestrain);
    elastic->SetPreStress(-44.3, -58.2, 0., -53.8);
    
    TPZGeoMesh *gmesh = new TPZGeoMesh;
    GeoMeshClass::WellBore2d(gmesh);
    ofstream arg("wellgeomeshlog.txt");
    gmesh->Print(arg);
    
    TPZTensor<STATE> Confinement;
    Confinement.XX() = -44.3;
    Confinement.YY() = -58.2;
    Confinement.ZZ() = -53.8;

    std::string vtkFile = "pocoelastic.vtk";
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(gmesh);
    compmesh1->InsertMaterialObject(elastic);
    CmeshWell(compmesh1,elastic,Confinement,-29.3);
    TPZAnalysis analysis(compmesh1,cout);
    TPZFStructMatrix structmatrix(compmesh1);
	analysis.SetStructuralMatrix(structmatrix);
    
	TPZSkylineStructMatrix full(compmesh1);
	analysis.SetStructuralMatrix(full);
    
	TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    analysis.Run();
    ////Post Processing
    TPZStack<std::string> vecnames,scalnames;

    
    vecnames.Push("displacement");
    vecnames.Push("Stress");
    scalnames.Push("SigmaX");
    scalnames.Push("SigmaY");
    scalnames.Push("SigmaZ");
    scalnames.Push("J2");
    scalnames.Push("I1");
    
    const int dim = 2;
    analysis.DefineGraphMesh(dim,scalnames,vecnames,vtkFile);
    analysis.PostProcess(0);
    

}



void RotationMatrix(TPZFMatrix<REAL> &R, double thetaRad, int axis)
{
    R.Resize(3,3);
    
    switch (axis)
    {
            
        case 0://ROTATE ABOUT X
            
            R.Put(0,0,1.);
            R.Put(1,1,cos(thetaRad));R.Put(1,2,sin(thetaRad));
            R.Put(2,1,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
            
            break;
            
        case 1://ROTATE ABOUT Y
            
            R.Put(1,1,1.);
            R.Put(0,0,cos(thetaRad));R.Put(0,2,sin(thetaRad));
            R.Put(2,0,-sin(thetaRad));R.Put(2,2,cos(thetaRad));
            
            break;
            
        case 2://ROTATE ABOUT Z
            
            R.Put(0,0,cos(thetaRad)); R.Put(0, 1,sin(thetaRad));
            R.Put(1,0,-sin(thetaRad)); R.Put(1, 1,cos(thetaRad));
            R.Put(2,2,1.);
            
            break;
            
        default:
            
            std::cout << " NON SPECIFIED AXIS "<<std::endl;
            break;
    }
    
}


void RotateMatrix(TPZFMatrix<REAL> &Mat, double thetaRad,int rotateaboutaxes)
{
    TPZFMatrix<REAL> R;
    TPZFMatrix<REAL> temp;
    RotationMatrix(R, thetaRad,rotateaboutaxes);
    
    if(R.Cols() != Mat.Rows())
    {
        cout << " \n -- MATRIX WITH INCOMPATIBLE SHAPES -- \n";
        DebugStop();
    }
    
    
    int matcol = Mat.Cols();
    TPZFMatrix<REAL> RT;
    if(matcol == 1)
    {
        R.Transpose(&RT);
        temp = RT*Mat;
        Mat = temp;
    }
    
    else
    {
        R.Transpose(&RT);
        temp = RT*Mat;
        Mat = temp*R;
    }
    
}

void RotateMesh(TPZGeoMesh &geomesh, REAL angleDegree,int rotateabout)
{
    REAL pi = M_PI;
    REAL th = angleDegree/180. * pi;
    
    for(int node = 0; node < geomesh.NodeVec().NElements(); node++)
    {
        TPZFMatrix<REAL> nodeCoord(3,1,0.);
        for(int c = 0; c < 3; c++)
        {
            nodeCoord.Put(c, 0, geomesh.NodeVec()[node].Coord(c));
        }
        RotateMatrix(nodeCoord,th,rotateabout);
        for(int c = 0; c < 3; c++)
        {
            geomesh.NodeVec()[node].SetCoord(c,nodeCoord(c,0));
        }
    }
}


void calcSDBar()
{
    REAL theta= 0;
    int axes=2;
    TPZGeoMesh *barmesh1;
    barmesh1 = BarMesh(0);
    RotateMesh(*barmesh1, theta, axes);
    ofstream arg("barmesh45.txt");
    barmesh1->Print(arg);
    
    TPZCompEl::SetgOrder(1);
	TPZCompMesh *compmesh1 = new TPZCompMesh(barmesh1);
	
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> SD;
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::McCormicRanchSand(SD);
    
	TPZMatElastoPlastic<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> > PlasticSD(1);
	PlasticSD.SetPlasticity(SD);
    
    TPZMaterial *plastic(&PlasticSD);
    compmesh1->InsertMaterialObject(plastic);
    
    Cmesh(compmesh1,plastic,theta,axes);
    
    ofstream arg2("barcmesh45.txt");
    compmesh1->Print(arg2);
    
	
	TPZElastoPlasticAnalysis EPAnalysis(compmesh1,cout);
    SolverSetUp2(EPAnalysis,compmesh1);
	
	TPZPostProcAnalysis PPAnalysis(compmesh1);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
    
    PPAnalysis.SetStructuralMatrix(structmatrix);
    TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariables2(PostProcVars,scalNames,vecNames);
    
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =3;
    
    
    PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"x.vtk");
    PPAnalysis.PostProcess(0/*pOrder*/);
	
    TPZFMatrix<REAL> BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix<REAL> val1(3,1,0.);TPZFMatrix<REAL> val2(3,1,0.);TPZFMatrix<REAL> BeginForce(3,1,0.);TPZFMatrix<REAL> EndForce(3,1,0.);
	
	int nsteps,taxa,nnewton;
	nnewton = 30;
	nsteps =10;
	taxa = 1;
	REAL beginforce = -0.0001;
	REAL endforce = -0.00145;
    int bc =-2;
    
    ManageIterativeProcess2(EPAnalysis ,beginforce,endforce,bc, nsteps,&PPAnalysis);
    EPAnalysis.Solution().Print("EPAnalysis.Solution");
    
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"x.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
    
    
    
}

void Cmesh(TPZCompMesh *CMesh, TPZMaterial * mat,REAL theta,int axes)
{
    
    TPZFMatrix<REAL> k1(3,3,0.);
    TPZFMatrix<REAL> f1(3,1,0.);
    RotateMatrix(k1,theta, axes);
    RotateMatrix(f1,theta, axes);
    TPZMaterial *bc1 = mat->CreateBC(mat,-1,0,k1,f1);
    CMesh->InsertMaterialObject(bc1);
    TPZFMatrix<REAL> k2(3,3,0.);
    TPZFMatrix<REAL> f2(3,1,0.);
    f2(0,0)=1.;
    RotateMatrix(k2,theta, axes);
    RotateMatrix(f2,theta, axes);
    TPZMaterial * bc2 = mat->CreateBC(mat,-2,1,k2,f2);
    CMesh->InsertMaterialObject(bc2);
    
    CMesh->AutoBuild();
    
}


TPZGeoMesh * BarMesh(int h)
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
    
    return gMesh;
    
}

void SolverSetUp2(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    
    //TPZFStructMatrix full(fCmesh)
	TPZSkylineStructMatrix full(fCmesh);
    full.SetNumThreads(0);
	an.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
    //  step.SetDirect(ELDLt);
    //  step.SetJacobi(5000, 1.e-12,0);
    step.SetDirect(ECholesky);
    // step.SetDirect(ELU);
	an.SetSolver(step);
	



}

void SetUPPostProcessVariables2(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	  scalnames.Resize(2);
	  scalnames[0] = "Alpha";
    scalnames[1] = "PlasticSqJ2";
//    scalnames[1] = "PlasticSteps";
//    scalnames[2] = "VolElasticStrain";
//    scalnames[3] = "VolPlasticStrain";
//    scalnames[4] = "VolTotalStrain";
//    scalnames[5] = "I1Stress";
//    scalnames[6] = "J2Stress";

    
	vecnames.Resize(2);
	vecnames[0] = "Displacement";
    vecnames[1] = "YieldSurface";
//    vecnames[1] = "NormalStress";
//    vecnames[2] = "ShearStress";
//    vecnames[3] = "NormalStrain";
//    vecnames[4] = "ShearStrain";
    
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

void ManageIterativeProcess2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost)
{
    
	TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    
	if(!pBC)return;
	
    
    REAL increment =  (valEnd-valBeg)/steps;
    cout<< "\n increment = " <<increment<<endl;
    TPZFMatrix<REAL> mattemp(3,1,0.);
    mattemp(0,0)=valBeg;
    //RotateMatrix(mattemp, 0,1);
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
        
        analysis.TransferSolution(*pPost);
        cout << "\nPostSolution "<< endl<< pPost->Solution();
        pPost->PostProcess(0);
        cout << pPost->Solution();
        
	}
}

void ManageIterativeProcessPesoProprio2(TPZElastoPlasticAnalysis &analysis ,REAL valBeg, REAL valEnd,int BCId, int steps,TPZPostProcAnalysis * pPost)
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
//
//#include "TPZYCModifiedMohrCoulomb.h"
//
//void MohrCoulombTestX()
//{
//
//
//
//
//	ofstream outfiletxt("DPInTriaxialCompression.txt");
//
//	TPZTensor<REAL> stress, strain, deltastress, deltastrain;
//
//    deltastress.XX() =-1.;
//    deltastress.XY() = 0.;
//    deltastress.XZ() = 0.;
//    deltastress.YY() = -1.;
//    deltastress.YZ() = 0.;
//    deltastress.ZZ() = -0.1;
//    stress = deltastress;
//
//    //3 2 1 cai no inner
//    //1 2 3 cai no inner
//    //2 3 1 cai no inner
//    //2 1 6 cai no outer
//    //2 1 0 cai MC
//
//    //	TPZFNMatrix<6*6> Dep(6,6,0.);
//    //    deltastrain.XX() = -0.0002;
//    //	deltastrain.XY() = 0.;
//    //	deltastrain.XZ() = 0.;
//    //	deltastrain.YY() = -0.00000001;
//    //	deltastrain.YZ() = 0.;
//    //	deltastrain.ZZ() = -0.00000001;
//    //	strain=deltastrain;
//
//    TPZPlasticStep<TPZYCModifiedMohrCoulomb,TPZThermoForceA,TPZElasticResponse> modifiedmohr;
//
//    REAL pi = M_PI;
//    REAL cohesion = 11.2033; //yield- coesao inicial correspondeno a fck igual 32 Mpa
//    REAL phi =  20./180. * pi; //phi=20
//    REAL hardening = 1000.; //Modulo de hardening da coesao equivante 1 Mpa a cada 0.1% de deformacao
//    REAL young = 20000.;
//    REAL poisson = 0.;
//    modifiedmohr.fYC.SetUp(phi);
//    modifiedmohr.fTFA.SetUp(cohesion, hardening);
//    modifiedmohr.fER.SetUp(young, poisson);
//
//	TPZMohrCoulomb Pstep;
//    Pstep.ConventionalConcrete(Pstep);
//
//    TPZDruckerPrager DP;
//    DP.ConventionalConcrete(DP,0);
//
//    int length =430;
//	for(int step=0;step<length;step++)
//	{
//		cout << "\nstep "<< step;
//
//		DP.ApplyLoad(stress,strain);
//        stress += deltastress;
//        //		cout<<  "\nstress " << stress << endl;
//        //		cout<<  "\nstrain " << strain << endl;
//        //		Pstep.Print(cout);
//        cout << "\n strain " << strain << endl;
//		TPZVec<REAL> phis(1);
//		DP.Phi(strain, phis);
//
//		cout << "\nphis " << phis << endl;
// 		//cout<<  "\nEigen " << eigenval << endl;
//        //if(step==50 || step == 100 ||step == 200 ||step == 300)deltastress*=-1.;
//		outfiletxt << fabs(strain.XX()) << " " << fabs(stress.XX()) << "\n";
//
//	}
//
//
//}
//




void SolverSet(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    
    //TPZFStructMatrix full(fCmesh)
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
    // step.SetDirect(ELDLt);
    //  step.SetJacobi(5000, 1.e-12,0);
    step.SetDirect(ECholesky);
    // step.SetDirect(ELU);
	an.SetSolver(step);
	
}

void PostProcessVariables(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
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


template <class T>
void WellboreLoadTest(stringstream & fileName, T & mat,
                      REAL loadMultipl, REAL plasticTol)
{
	REAL L, LDA, alpha, GammaSea, MudWeight, kh, kH, OCR, s;
	s = loadMultipl;
	REAL SigmaV=0., SigmaH=0., Sigmah=0.;
	REAL Pa = 14.7;
	int ncirc, ioRatio, pOrder, valType;
	
	cout << "\nMesh data: ncirc?(int) ";
	cin >> ncirc;
	
	cout << "Mesh data: ioratio?(sugg. 10.) ";
	cin >> ioRatio;
	
	cout << "Mesh data: pOrder? ";
	cin >> pOrder;
	
	fileName << "_ncirc" << ncirc << "_IO" << ioRatio << "_p" << pOrder;
	
	cout << "\nSelect one of the following load case (Load):";
	cout << "\n0) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.8 kH=0.9 OCR=1.10";
	cout << "\n1) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.6 kH=0.8 OCR=1.10";
	cout << "\n2) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=0.5";
	cout << "\n3) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=1.0";
	cout << "\n4) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=9.5 kh=0.80 kH=1.30 OCR=1.10 Biot=1.0";
	cout << "\n";
	
	cin >> valType;
    
	switch(valType)
	{
		case(0):
		    L=3000;
		    LDA=1200;
		    alpha = 0.8;
		    GammaSea = 8.6;
		    MudWeight = 9.2;
		    kh = 0.8;
		    kH = 0.9;
		    OCR = 1.1;
		    fileName << "_Load0";
            break;
		case(1):
		    L=3000;
		    LDA=1200;
		    alpha = 0.8;
		    GammaSea = 8.6;
		    MudWeight = 9.2;
		    kh = 0.6;
		    kH = 0.8;
		    OCR = 1.1;
		    fileName << "_Load1";
            break;
		case(2):
		    L=5227;
		    LDA=2135;
		    alpha = 0.5;
		    GammaSea = 9.5;
		    MudWeight = 9.5;
		    SigmaV= loadMultipl * 84.4 * 145.03773801/Pa;
		    kh = 0.96;
		    kH = 1.09;
		    OCR = 1.1;
		    fileName << "_Load2";
            break;
		case(3):
		    L=5227;
		    LDA=2135;
		    alpha = 1.;
		    GammaSea = 9.5;
		    MudWeight = 10;
		    SigmaV= loadMultipl * 84.4 * 145.03773801/Pa;
		    kh = 0.96;
		    kH = 1.09;
		    OCR = 1.1;
		    fileName << "_Load3";
            break;
		case(4):
		    L=5227;
		    LDA=2135;
		    alpha = 1.;
		    GammaSea = 9.5;
		    MudWeight = 10;
		    SigmaV= loadMultipl * 84.4 * 145.03773801/Pa;
		    kh = 0.80;
		    kH = 1.30;
		    OCR = 1.1;
		    fileName << "_Load3";
            break;
		default:
		    L=0;
		    LDA=0;
		    alpha = 0.;
		    GammaSea = 0.;
		    MudWeight = 0.;
		    kh = 0.;
		    kH = 0.;
		    OCR = 0.;
			cout << "Unhandled Case. Exiting...";
            break;
	}
	
	fileName << ".vtk";
	
	cout << endl << fileName.str() << endl;
	
	//Overconsolidating the material
	
	TPZTensor<REAL> OCStress, beginOCStress, loadStress, loadStress2, initialStrain, FarFieldStress, TestStress;
	TPZFNMatrix<3*3> BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFNMatrix<3*1> val1(3,1,0.);
	
	const REAL a = 0.17046;
	REAL PorePressure = s * GammaSea * L * a / Pa;
    
    cout << PorePressure;
    cout <<"<\n>"<<SigmaV;
	if(SigmaV==0.) SigmaV  = s * (0.9 * (L-LDA) / 0.3048 + GammaSea * LDA * a) / Pa;
	Sigmah       = kh * (SigmaV - PorePressure * alpha) + PorePressure * alpha;
	SigmaH       = kH * (SigmaV - PorePressure * alpha) + PorePressure * alpha;
    
	REAL FluidWeight  = s * MudWeight * a * L / Pa;
	
	FarFieldStress.fData[_XX_] = Sigmah - PorePressure * alpha;
	FarFieldStress.fData[_YY_] = SigmaH - PorePressure * alpha;
	FarFieldStress.fData[_ZZ_] = SigmaV - PorePressure * alpha;
	FarFieldStress *= Pa;
    
	beginOCStress.Identity();
	beginOCStress *= s * 0.01 *Pa;
	
	OCStress = FarFieldStress;
	OCStress *= OCR;
	
	loadStress.Identity();
	
	loadStress *= (Sigmah - PorePressure * alpha) *Pa;
	
	loadStress2.Identity();
	loadStress2 *= (FluidWeight - PorePressure * alpha) *Pa;
	
	loadStress.    CopyToTensor(EndStress);
	loadStress2.   CopyToTensor(EndStress2);
	FarFieldStress.CopyToTensor(BeginStress);
	
	cout << "\nInitial Stress State: " << FarFieldStress;
	cout << "\nLoad Stress State  : " << loadStress;
	cout << "\nLoad Stress State 2: " << loadStress2;
	cout << "\n" ;
    
	PrepareInitialMat(mat, beginOCStress, OCStress, 10);
	
	// Returning the strain state back to the correspondent imposed stress state
    
	mat.ApplyLoad(FarFieldStress, initialStrain);
	
	cout << "\nApplied Desired Stress State: " << FarFieldStress <<"\n resulted in strain: "<< initialStrain <<"\n";
	
	mat.ApplyStrainComputeSigma(initialStrain, TestStress);
	
	cout << "\nApplied Desired Strain State: " << initialStrain <<"\n resulted in stress: "<< TestStress <<"\n";
    
	cout << "\n Plastic State = " << mat.GetState();
	
	//Attributing the material history to the PZ ElastoPlastic material object
	
	TPZMatElastoPlastic<T> EPMat(1);
	
	EPMat.SetPlasticity(mat);
	
    TPZCompMesh * pCMesh;
    
    
    //PZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh);
    //TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
    
    pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
	
    //End of material initialization
	//TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
	
	//TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
	
	//building analysis
	TPZElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
    
    
    SolverSet(EPAnalysis,pCMesh);
    
    
	TPZPostProcAnalysis PPAnalysis(pCMesh);
    
    
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
    
    
    
	
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
    
	PostProcessVariables(PostProcVars,scalNames, vecNames);
    
    PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	cout << "\nTransfering initial Solutions\n";
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	
	PPAnalysis.DefineGraphMesh(3,scalNames,vecNames,fileName.str());
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	
	EPAnalysis.ManageIterativeProcess(cout, 1.e-5, 10,
									  -7 /*BCId*/, 2 /*nsteps*/, 1/*PGRatio*/,
									  BeginStress/*val1Begin*/, EndStress/*val1End*/,
									  val1/*val2Begin*/, val1/*val2End*/,
									  &PPAnalysis, pOrder);
	
	EPAnalysis.ManageIterativeProcess(cout, 1.e-5, 10,
									  -7 /*BCId*/, 2 /*nsteps*/, 1/*PGRatio*/,
									  EndStress/*val1Begin*/, EndStress2/*val1End*/,
									  val1/*val2Begin*/, val1/*val2End*/,
									  &PPAnalysis, pOrder);
	
	
	cout << "\nProblem Solved. Accepting new Solutions\n";
	
	cout << "\nClosing Mesh\n";
	
	PPAnalysis.CloseGraphMesh();
	
	cout << "\nExiting\n";
    
    
	
    return;
    
    
}



template <class T>
void PorousWellboreLoadTest(stringstream & fileName, T & mat,
                            REAL loadMultipl, REAL plasticTol)
{
	REAL L, LDA, alpha, GammaSea, MudWeight, kh, kH, OCR, s;
	REAL perm, mu, storageEps, rhof, SigmaV=0., SigmaH=0., Sigmah=0.;
	REAL Pa = 14.7;
	s = loadMultipl;
	int ncirc, ioRatio, pOrder, valType;
	
	cout << "\nMesh data: ncirc?(int) ";
	cin >> ncirc;
	
	cout << "Mesh data: ioratio?(sugg. 10.) ";
	cin >> ioRatio;
	
	cout << "Mesh data: pOrder? ";
	cin >> pOrder;
	
	fileName << "_ncirc" << ncirc << "_IO" << ioRatio << "_p" << pOrder;
	
	cout << "\nSelect one of the following load case (Load):";
	cout << "\n0) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.8 kH=0.9 OCR=1.10";
	cout << "\n1) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.6 kH=0.8 OCR=1.10";
	cout << "\n2) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=0.5";
	cout << "\n3) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=1.0";
	cout << "\n";
	
	cin >> valType;
    
	switch(valType)
	{
		case(0):
		    L=3000;
		    LDA=1200;
		    alpha = 0.8;
		    GammaSea = 8.6;
		    MudWeight = 9.2;
		    kh = 0.8;
		    kH = 0.9;
		    OCR = 1.1;
		    fileName << "_Load0";
            break;
		case(1):
		    L=3000;
		    LDA=1200;
		    alpha = 0.8;
		    GammaSea = 8.6;
		    MudWeight = 9.2;
		    kh = 0.6;
		    kH = 0.8;
		    OCR = 1.1;
		    fileName << "_Load1";
            break;
		case(2):
		    L=5227;
		    LDA=2135;
		    alpha = 0.5;
		    GammaSea = 9.5;
		    MudWeight = 10;
		    SigmaV= loadMultipl * 84.4 * 145.03773801/Pa;
		    kh = 0.96;
		    kH = 1.09;
		    OCR = 1.1;
		    fileName << "_Load2";
            break;
		case(3):
		    L=5227;
		    LDA=2135;
		    alpha = 1.;
		    GammaSea = 9.5;
		    MudWeight = 10;
		    SigmaV= loadMultipl * 84.4 * 145.03773801/Pa;
		    kh = 0.96;
		    kH = 1.09;
		    OCR = 1.1;
		    fileName << "_Load3";
            break;
		default:
		    L=0;
		    LDA=0;
		    alpha = 0.;
		    GammaSea = 0.;
		    MudWeight = 0.;
		    kh = 0.;
		    kH = 0.;
		    OCR = 0.;
			cout << "Unhandled Case. Exiting...";
            break;
	}
	
	perm = 2.961e-13; //[m2]
	mu = 2.e-9; // [MPa.s]
	storageEps = 1.156e-4;// [MPa^-1]//0.00000000001;
	rhof = 0.; // null fluid density in order to ignore gravity effects.
	
	fileName << ".vtk";
	
	cout << endl << fileName.str() << endl;
	
	//Overconsolidating the material
	
	TPZTensor<REAL> OCStress, beginOCStress, loadStress, loadStress2,
    initialStrain, FarFieldStress, EffFarFieldStress, TestStress;
	TPZFNMatrix<3*3> BeginStress(4,4,0.), EndStress(4,4,0.), EndStress2(4,4,0.);
	TPZFNMatrix<3*1> val1(4,1,0.);
	
	const REAL a = 0.17046;
	REAL PorePressure = s * GammaSea * L * a / Pa;
	if(SigmaV==0.)SigmaV  = s * (0.9 * (L-LDA) / 0.3048 + GammaSea * LDA * a) / Pa;
	Sigmah       = kh * (SigmaV - PorePressure * alpha) + PorePressure * alpha;
	SigmaH       = kH * (SigmaV - PorePressure * alpha) + PorePressure * alpha;
    
	REAL FluidWeight  = s * MudWeight * a * L / Pa;
	
	FarFieldStress.fData[_XX_] = Sigmah;
	FarFieldStress.fData[_YY_] = SigmaH;
	FarFieldStress.fData[_ZZ_] = SigmaV;
	
	EffFarFieldStress.fData[_XX_] = Sigmah - PorePressure * alpha;
	EffFarFieldStress.fData[_YY_] = SigmaH - PorePressure * alpha;
	EffFarFieldStress.fData[_ZZ_] = SigmaV - PorePressure * alpha;
	FarFieldStress *= Pa;
	EffFarFieldStress *= Pa;
    
	beginOCStress.Identity();
	beginOCStress *= s * 0.01 *Pa;
	
	OCStress = EffFarFieldStress;
	OCStress *= OCR;
	
	loadStress.Identity();
	
	loadStress *= Sigmah * Pa;
	
	loadStress2.Identity();
	loadStress2 *= FluidWeight * Pa;
	
	loadStress.    CopyToTensor(EndStress);
	EndStress.Resize(4,4); // reserving space for the PorePressure BC
	loadStress2.   CopyToTensor(EndStress2);
	EndStress2.Resize(4,4);
	FarFieldStress.CopyToTensor(BeginStress);
	BeginStress.Resize(4,4);
	
	cout << "\nInitial Total Stress State: " << FarFieldStress;
	cout << "\nLoad Total Stress State  : " << loadStress;
	cout << "\nLoad Total Stress State 2: " << loadStress2;
	cout << "\n" ;
    
	PrepareInitialMat(mat, beginOCStress, OCStress, 10);
	/*
     TPZPlasticState<REAL> state;
     state.fEpsT.fData[0] = -0.0212118;
     state.fEpsT.fData[3] = -0.0230264;
     state.fEpsT.fData[5] = -0.024841;
     state.fEpsP.fData[0] = -0.0189783;
     state.fEpsP.fData[3] = -0.0204523;
     state.fEpsP.fData[5] = -0.0219263;
     state.fAlpha = -0.0613569;
     
     mat.SetState(state);
     */
	// Returning the strain state back to the correspondent imposed stress state
    
	mat.ApplyLoad(EffFarFieldStress, initialStrain);
	
	cout << "\nApplied Desired Stress State: " << EffFarFieldStress <<"\n resulted in strain: "<< initialStrain <<"\n";
	
	mat.ApplyStrainComputeSigma(initialStrain, TestStress);
	
	cout << "\nApplied Desired Strain State: " << initialStrain <<"\n resulted in stress: "<< TestStress <<"\n";
    
	cout << "\n Plastic State = " << mat.GetState() << "\n";
	
	//Attributing the material history to the PZ ElastoPlastic material object
	
	TPZMatPorous<T> EPMat(1);
	
	EPMat.SetUp(perm, mu, storageEps, alpha, rhof);
	EPMat.SetPlasticity(mat);
	EPMat.SetPorePressure(fabs(PorePressure * Pa));
	
    //TPZPoroElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
	
	
    TPZCompMesh * pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
	
    
    //End of material initialization
	//TPZPoroElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
	
    
	
	{// Replacing BC to include Pore Pressure
		TPZMaterial * mat(&EPMat);
		TPZFMatrix<REAL> val2(4,1,0.), val1(4,4,0.);
		val2.Zero();
		val2(0,0) = 1.;
		val2(1,0) = 1.;
		val2(3,0) = 1.;
		val1(3,3) = 0;// no pore pressure change
		TPZMaterial * bc;
		bc = mat->CreateBC(mat, -5, 3, val1, val2);	//Directional Dirichlet BCType
		pCMesh->InsertMaterialObject(bc);
	}
	
	//building analysis
	TPZPoroElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
	
	//EPAnalysis.SetBiCGStab(50000, 1.e-10);
    SolverSet(EPAnalysis, pCMesh);
    
    
    
    
    
	// Preparing Post Process
    
	TPZPostProcAnalysis PPAnalysis(pCMesh);
    
    TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
	
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	PostProcessVariables(PostProcVars,scalNames, vecNames);
    
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	cout << "\nTransfering initial Solutions\n";
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	
	PPAnalysis.DefineGraphMesh(3,scalNames,vecNames,fileName.str());
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
	cout << "\nInitial Solution Exported. Solving Problem\n";
	
	EPAnalysis.SetDeltaT(1.e-12);
	
	EPAnalysis.ManageIterativeProcess(cout, 1.e-4, 10,
									  -7 /*BCId*/, 5 /*nsteps*/, 1/*PGRatio*/,
									  BeginStress/*val1Begin*/, EndStress/*val1End*/,
									  val1/*val2Begin*/, val1/*val2End*/,
									  &PPAnalysis, pOrder);
	
	EPAnalysis.ManageIterativeProcess(cout, 1.e-4, 10,
									  -7 /*BCId*/, 5 /*nsteps*/, 1/*PGRatio*/,
									  EndStress/*val1Begin*/, EndStress2/*val1End*/,
									  val1/*val2Begin*/, val1/*val2End*/,
									  &PPAnalysis, pOrder);
	REAL time = 0, deltaT = 0;
	for(int i = 0; i < 10; i++)
	{
		REAL currentTime = pow(10.,i);
		deltaT = currentTime - time;
		time = currentTime;
		cout << "\n Evoluting to time " << currentTime;
		
		EPAnalysis.SetDeltaT(deltaT);
        
		EPAnalysis.Run(cout, 1.e-4, 10, &PPAnalysis, pOrder);
	}
	
	//PPAnalysis.PostProcess(pOrder);
	
	cout << "\nProblem Solved. Accepting new Solutions\n";
	
	cout << "\nClosing Mesh\n";
	
	PPAnalysis.CloseGraphMesh();
	
	cout << "\nExiting\n";
    
	
    return;
    
    
}

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

int main ()
{
  
//int startfrom = 0;
startfrom = 4;
    
#ifdef MACOS
    
    ENABLE_FPO_EXCEPTIONS;
    ATTACH_FPO_SIGNAL;
    
#endif
    
    InitializePZLOG();
    gRefDBase.InitializeAllUniformRefPatterns();
	
   // TPZGeoMesh * mesh = GeoMeshClass::WellBore2d();
    
    /* QuarterWellboreGeom2d(int ncirc,
     REAL ioratio,
     TPZVec< TPZVec<REAL> > &pt,
     TPZVec< TPZVec<int> > &el,
     TPZVec< MElementType > &elType,
     TPZVec< TPZString > &elName,
     int & nrad)
     
     */
//    TPZWellBoreAnalysis::CheckDeformation();
//    int startfrom = 3;
    TPZWellBoreAnalysis well;
    REAL innerradius = 4.25*0.0254;
    REAL outerradius = 3.;
    REAL computedquarter = 7.05761678496926; 
    REAL sqj2_refine = 0.0007;
    std::cout << std::setprecision(15);
    if (startfrom == 0) 
    {
        well.SetInnerOuterRadius(innerradius, outerradius);
        TPZManVector<STATE,3> confinement(3,0.);
        confinement[0] = -45.9;
        confinement[1] = -62.1;
        confinement[2] = -48.2;
//        well.SetConfinementStresses(confinement, 28.9);
//
        REAL effectivePressure = 19.5; // 19.5 ou 23.4 ou 28.9 
        well.SetConfinementStresses(confinement, effectivePressure);
        REAL poisson = 0.203;
        REAL elast = 29269.;
        REAL A = 152.54;
        REAL B = 0.0015489;
        REAL C = 146.29;
        REAL R = 0.91969;
        REAL D = 0.018768;
        REAL W = 0.006605;
        well.SetSanderDiMaggioParameters(poisson, elast, A, B, C, R, D, W);
        int divisions = 20;
        REAL delx = 0.2*innerradius*M_PI_2/divisions;
        TPZManVector<int,2> numdiv(2,divisions);
        numdiv[1] = 40;
        well.SetMeshTopology(delx, numdiv);
        well.GetCurrentConfig()->CreateMesh();
        int porder = 2;
        well.GetCurrentConfig()->CreateComputationalMesh(porder);
        well.GetCurrentConfig()->CreatePostProcessingMesh();
//        well.LinearConfiguration(1);
        well.PostProcess(0);
        TPZBFileStream save;
        save.OpenWrite("Wellbore0.bin");
        well.Write(save);
    }    
//    if (startfrom == 0) 
//    {
//        TPZWellBoreAnalysis::StandardConfiguration(well);
//        TPZBFileStream save;
//        save.OpenWrite("Wellbore0.bin");
//        well.Write(save);
//    }
    if (startfrom ==1)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore0.bin");
        well.Read(read);
    }
    if (startfrom <= 1)
    {
        int nsteps = 3;
        int numnewton = 80;
        well.ExecuteInitialSimulation(nsteps, numnewton);
        TPZBFileStream save;
        save.OpenWrite("Wellbore1.bin");
        well.Write(save);

    
    }

    if (startfrom == 2)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore1.bin");
        well.Read(read);
    }
    if (startfrom <= 2) {
        
        
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        //        well.PostProcess(1);
//        int nsteps = 3;
//        int numnewton = 80;
        //        well.ExecuteSimulation();
        //        well.ExecuteInitialSimulation(nsteps,numnewton);
        
        //well.PRefineElementAbove(0., 2);
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation();
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.00000001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "Analytical - computed area " << analyticarea-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        std::cout << "Saving Wellbore2.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore2.bin");
        well.Write(save);
    }
    if (startfrom == 3) {
        TPZBFileStream read;
        read.OpenRead("Wellbore2.bin");
        well.Read(read);        
    }
    if (startfrom <= 3)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.010;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.79;
        //vvalor de a e b para sqJ2 = 0.0005 p 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.03409;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.829545;
        //vvalor de a e b para sqJ2 = 0.0007 p 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.014;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.90;
        std::multimap<REAL, REAL> polygonalChain;
        well.GetJ2Isoline(0.0007, polygonalChain);
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),coefs(2),verify(2);
        ellips.Getcoefficients(coefs);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        verify[0] = ratios[0]/well.GetCurrentConfig()->fInnerRadius;
        verify[1] = ratios[1]/well.GetCurrentConfig()->fInnerRadius;
        REAL a = ratios[0];
        REAL b = ratios[1];
        well.AddEllipticBreakout(a, b);
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        //well.PostProcess(1);
        well.ExecuteSimulation();
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation();
        
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.0001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        

        std::cout << "Saving Wellbore3.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore3.bin");
        well.Write(save);
        
        if(startfrom == 3)
        {
            TPZFileStream write;
            write.OpenWrite("afterread.txt");
            well.GetCurrentConfig()->Write(write);
        }
        else {
            TPZFileStream write;
            write.OpenWrite("afterrun.txt");
            well.GetCurrentConfig()->Write(write);
        }
        
    }
    if (startfrom == 4)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore3.bin");
        well.Read(read);        
    }
    if (startfrom <= 4)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.040;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.46;
        //vvalor de a e b para sqJ2 = 0.0005 p 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.10227;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.596591;
        //vvalor de a e b para sqJ2 = 0.0007 Pef = 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.044;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.66;
        std::multimap<REAL, REAL> polygonalChain;
        well.GetJ2Isoline(0.0007, polygonalChain);
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),coefs(2),verify(2);
        ellips.Getcoefficients(coefs);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        well.AddEllipticBreakout(ratios[0], ratios[1]);
        well.PostProcess(1);
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.ExecuteSimulation();
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation();
        
        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.000001);
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        

        std::cout << "Saving Wellbore4.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore4.bin");
        well.Write(save);
    }
    if (startfrom == 5) {
        TPZBFileStream read;
        read.OpenRead("Wellbore4.bin");
        well.Read(read);        
    }
    if (startfrom <= 5)
    {
        // valor de a e b para sqJ2 = 0.00025 E USANDO Pef = 23.4
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.082;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.30;
        // valor de a e b para sqJ2 = 0.0005 P = 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.17045;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.465909;
        // valor de a e b para sqJ2 = 0.0007 P = 19.5
//        REAL a = well.GetCurrentConfig()->fInnerRadius*1.099;
//        REAL b = well.GetCurrentConfig()->fInnerRadius*0.419;
        std::multimap<REAL, REAL> polygonalChain;
        well.GetJ2Isoline(0.0007, polygonalChain);
        TPZProjectEllipse ellips(polygonalChain);
        TPZManVector<REAL,2> center(2),ratios(2),coefs(2),verify(2);
        ellips.Getcoefficients(coefs);
        ellips.StandardFormatForSimpleEllipse(center, ratios);
        well.AddEllipticBreakout(ratios[0], ratios[1]);
        well.GetCurrentConfig()->ModifyWellElementsToQuadratic();
        well.PostProcess(1);
        well.ExecuteSimulation();
        well.DivideElementsAbove(0.0001);
        well.PRefineElementAbove(0.0001, 3);
        well.ExecuteSimulation();

        REAL analyticarea = M_PI*(outerradius*outerradius-innerradius*innerradius)/4.;
        REAL originalarea = well.GetCurrentConfig()->ComputeTotalArea();
        REAL plastifiedarea = well.GetCurrentConfig()->ComputeAreaAboveSqJ2(sqj2_refine);
        REAL openingangle = well.GetCurrentConfig()->OpeningAngle(0.0001);
        std::cout << "Analytical area " << analyticarea << " computed area " << originalarea << std::endl;
        std::cout << "computed quarter - domain area " << computedquarter-originalarea << std::endl;
        std::cout << "Plastified area " << plastifiedarea << std::endl;
        std::cout << "Opening angle " << openingangle << std::endl;
        
        
        std::cout << "Saving Wellbore5.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore5.bin");
        well.Write(save);
    }
    if (startfrom == 6)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore5.bin");
        well.Read(read);        
        
    }
    return 0;
    if (startfrom <= 6)
    {
//        well.VerifyGlobalEquilibrium();
        if(0)
        {
            
            TPZStack<std::string> postprocess;
            postprocess.Push("I1J2Stress");
            TPZFMatrix<STATE> valuetable;
            //TPZManVector<REAL,3> x(3,0.);
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(well.GetCurrentConfig()->fCMesh.ElementVec()[0]);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &intpoints = intel->GetIntegrationRule();
            TPZManVector<REAL> ksi(2,0.),xco(3,0.);
            REAL weight;
            intpoints.Point(0, ksi, weight);
            TPZMaterialData data;
            intel->InitMaterialData(data);
            intel->ComputeRequiredData(data, ksi);
            TPZManVector<long> memindices(intpoints.NPoints());
            intel->GetMemoryIndices(memindices);
            data.intGlobPtIndex = memindices[0];
            TPZMaterial *mat = intel->Material();
            int varindex = mat->VariableIndex("I1J2Stress");
            int nvar = mat->NSolutionVariables(varindex);
            TPZManVector<STATE> post(nvar);
            mat->Solution(data, varindex, post);
            std::cout << "Post processed " << post << std::endl;
            intel->Reference()->X(ksi, xco);
            //x[0] = 1.1;
            well.PostProcessedValues(xco , postprocess, valuetable);
            valuetable.Print("Post processed I1=J2",std::cout);
            xco.Fill(0.);
            well.PostProcessedValues(xco , postprocess, valuetable);
            valuetable.Print("Post processed I1=J2",std::cout);
        }
        
         //well.ChangeMaterialId(-2, -6);
        well.DeleteElementsAbove(0.0005);
        well.ChangeMaterialId(-6, -2);
        well.ExecuteSimulation();
        std::cout << "Saving Wellbore7.bin\n";
        TPZBFileStream save;
        save.OpenWrite("Wellbore7.bin");
        well.Write(save);

    }
    if (startfrom == 7)
    {
        TPZBFileStream read;
        read.OpenRead("Wellbore7.bin");
        well.Read(read);
    }
    return 0;
    
    
	int testNumber, matNumber;
	TPZPlasticBase *pMat;
	TPZLadeKim * pLK = NULL;
	TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> * pSD = NULL;
	typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
	TPZDruckerPrager * pDP = NULL;
	REAL loadMultipl;
	stringstream fileName, copyStr;
	REAL plasticTol;
	
	cout << "\nPlease enter test type:";
	cout << "\n0) Wellbore Drilling Load";
	cout << "\n1) Wellbore Drilling Load - Porous Medium";
	cout << "\n2) Graphical representation of Yield surface";
    
	//cin >> testNumber;
    
	
	cout << "\nMaterial Type:";
	cout << "\n1)Sandler Dimaggio: McCormic Ranch Sand";
	cout << "\n2)Sandler Dimaggio: McCormic Ranch Sand Mod";
	cout << "\n3)Sandler Dimaggio: Unconsolidated Deep Reservoir Sandstone [psi]";
	cout << "\n4)Sandler Dimaggio: Unconsolidated Deep Reservoir Sandstone [MPa]";
	cout << "\n5)Sandler Dimaggio: PRSMat [MPa]";
	cout << "\n6)Drucker Prager (Inscr MC): PRSMat [MPa]";
	cout << "\n7)Drucker Prager (Circunscr MC): PRSMat [MPa]";
	cout << "\n";
    
    //    cin >> matNumber;
    matNumber = 5;
	
	switch(matNumber)
	{
        case(1):
            pSD = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>();
		    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::McCormicRanchSand(*pSD);
			pMat = pSD;
		    fileName << "_SDMc";
		    loadMultipl = -0.001;
            break;
		case(2):
			pSD = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>();
		    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::McCormicRanchSandMod(*pSD);
			pMat = pSD;
		    fileName << "_SDMM";
		    loadMultipl = -0.001;
            break;
		case(3):
			pSD = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>();
		    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::UncDeepSandResPSI(*pSD);
			pMat = pSD;
		    fileName << "_SDDS";
		    loadMultipl = -1;
            break;
		case(4):
			pSD = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>();
		    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::UncDeepSandResMPa(*pSD);
			pMat = pSD;
		    fileName << "_SDDSMPa";
		    loadMultipl = -1/145.03773801;
            break;
		case(5):
			pSD = new TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>();
		    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1>::PRSMatMPa(*pSD);
			pMat = pSD;
		    fileName << "_PRSLMPa";
		    loadMultipl = -1/145.03773801;
            break;
		case(6):
            pDP = new TPZDruckerPrager();
			pDP->fYC.SetUp(/*phi*/ 29.7/180. * pi ,/*innerMCFit*/1);
			pDP->fTFA.SetUp(/*yield- coesao inicial*/ 12.8, /*k Modulo de hardening da coesao equivante 10^-3 Mpa a cada 0.1% de deformacao */1.);
			pDP->fER.SetUp(/*young*/ 29269., /*poisson*/ 0.203);
			pMat = pDP;
		    fileName << "_PRDPInscMPa";
		    loadMultipl = -1/145.03773801;
            break;
		case(7):
            pDP = new TPZDruckerPrager();
			pDP->fYC.SetUp(/*phi*/ 29.7/180. * pi ,/*innerMCFit*/0);
			pDP->fTFA.SetUp(/*yield- coesao inicial*/ 12.8, /*k Modulo de hardening da coesao equivante 10^-3 Mpa a cada 0.1% de deformacao */1.);
			pDP->fER.SetUp(/*young*/ 29269., /*poisson*/ 0.203);
			pMat = pDP;
		    fileName << "_PRDPCircMPa";
		    loadMultipl = -1/145.03773801;
            break;
		default:
			cout << "\nUnhandled Material Type. Exiting...";
            break;
	}
	
	cout << "\nPlastic Integration Tolerance:(sugg. 0.0001) ";
	
    //	cin >> plasticTol;
    plasticTol = 0.0001;
    
	fileName << "_pTol" << plasticTol;
    
    testNumber = 1;
    InitializePZLOG();
	
	switch(testNumber)
	{
		case(0):
			copyStr << fileName.str();
			fileName.str(""); // clearing the buffer
		    fileName << "WB" << copyStr.str();
		    if(pLK)WellboreLoadTest(fileName, *pLK, loadMultipl, plasticTol);
		    if(pSD)WellboreLoadTest(fileName, *pSD, loadMultipl, plasticTol);
		    if(pDP)WellboreLoadTest(fileName, *pDP, loadMultipl, plasticTol);
            break;
		case(1):
			copyStr << fileName.str();
			fileName.str(""); // clearing the buffer
		    fileName << "PorousWB" << copyStr.str();
		    if(pLK)PorousWellboreLoadTest(fileName, *pLK, loadMultipl, plasticTol);
		    if(pSD)PorousWellboreLoadTest(fileName, *pSD, loadMultipl, plasticTol);
		    if(pDP)PorousWellboreLoadTest(fileName, *pDP, loadMultipl, plasticTol);
            break;
        case 2:
        {
            
            
        }
        case 3://FIGURA 11-a
        {
            TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler;
            sandler.McCormicRanchSand(sandler);
            ofstream outfiletxty("FIGURA11A.txt");
            TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
            
            sandler.SetIntegrTol(0.00001);
            
            deltaeps.XX()= -0.0013;
            deltaeps.YY()=0;
            deltaeps.ZZ()=0;
            eps=deltaeps;
            
            for(int i=0;i<100;i++)
            {
                
                sandler.ApplyStrainComputeSigma(eps, sigma);//UCS
                if(i==49)
                {
                    deltaeps*=-1;
                }
                
                REAL sqrJ2 = sqrt(sigma.J2());
                REAL I1=sigma.I1();
                outfiletxty << -I1<< " " << sqrJ2 << "\n";
                eps+=deltaeps;
            }
            
            
        }
        case 4://FIGURA 12
        {
            TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler;
            sandler.McCormicRanchSand(sandler);
            ofstream outfiletxty("FIGURA12.txt");
            TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
            TPZPlasticState<REAL> state;
            
            sandler.SetIntegrTol(1);
            
            deltaeps.XX()= -0.0013;
            deltaeps.YY()=0;
            deltaeps.ZZ()=0;
            eps=deltaeps;
            
            for(int i=0;i<100;i++)
            {
                
                sandler.ApplyStrainComputeSigma(eps, sigma);//UCS
                cout << "\n sigma "<<sigma <<endl;
                REAL j2 = 0.;
                j2 = sigma.J2();
                state  = sandler.GetState();
                if(i==50)
                {
                    deltaeps*=-1;
                }
                
                outfiletxty << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
                
                eps+=deltaeps;
            }
            
            
        }
        case 5://FIGURA 13
        {
            TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler3;
            ofstream outfiletxt("FIGURA13.txt");
            TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
            
            deltasigma.XX()=-0.004;
            deltasigma.YY()=deltasigma.XX()*0.4;
            deltasigma.ZZ()=deltasigma.YY();
            sigma=deltasigma;
            
            sandler3.McCormicRanchSand(sandler3);
            sandler3.SetIntegrTol(0.0001);
            for(int i=0;i<100;i++)
            {
                sandler3.ApplyLoad(sigma,eps);
                outfiletxt << fabs(eps.XX()) << " " << fabs(sigma.XX()) << "\n";
                //  outfiletxtS << sigma.I1() << " " << sqrt(sigma.J2()) << "\n";
                if(i==100 || i==200)
                {
                    deltasigma*=-1;
                }
                sigma+=deltasigma;
            }
            
        }
            
        case 6:
        {
            
            TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> sandler3;
            ofstream outfiletxt("FIGURA14.txt");
            TPZTensor<REAL> deltaeps,eps,sigma,deltasigma;
            
            deltasigma.XX()=-0.004;
            deltasigma.YY()=deltasigma.XX()*0.8;
            deltasigma.ZZ()=deltasigma.YY();
            sigma=deltasigma;
            
            sandler3.McCormicRanchSand(sandler3);
            sandler3.SetIntegrTol(0.0001);
            for(int i=0;i<25;i++)
            {
                sandler3.ApplyLoad(sigma,eps);
                outfiletxt << fabs(sigma.XX()-sigma.YY()) << " " << fabs(eps.XX()) << "\n";
                sigma+=deltasigma;
            }
            
            
        }
            
            break;
		default:
			cout << "\nUnhandled Test Type. Exiting...";
            delete pMat;
	}
	
    return EXIT_SUCCESS;
    
    
}

#include "TPZGenSpecialGrid.h"

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD);

void VisualizeSandlerDimaggio(std::stringstream &fileName, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD)
{
    TPZVec<REAL> coords(0);
    TPZGeoMesh *gmesh = TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(coords, 0.001,1);
    TPZCompMesh *cgrid = new TPZCompMesh(gmesh);
    TPZManVector<STATE> force(3,0.);
    
    for (int imat=1; imat<4; imat++) {
        TPZMat2dLin *mat = new TPZMat2dLin(imat);
        TPZFNMatrix<9,STATE> xk(3,3,0.),xc(3,3,0.),xf(3,1,0.);
        mat->SetMaterial(xk,xc,xf);
        //    TPZMaterial * mat = new TPZElasticity3D(1,1.e5,0.2,force);
        cgrid->InsertMaterialObject(mat);
    }
    cgrid->AutoBuild();
    TPZFMatrix<REAL> elsol(cgrid->NElements(),1,0.);
    
    cgrid->ElementSolution() = elsol;
    
    TPZAnalysis an(cgrid);
    std::stringstream vtkfilename;
    vtkfilename << fileName.str();
    vtkfilename << ".vtk";
    std::ofstream meshout(vtkfilename.str().c_str());
	TPZVTKGeoMesh::PrintGMeshVTK(gmesh,meshout);
    BuildPlasticSurface(cgrid,pSD);
    TPZStack<std::string> scalnames,vecnames;
    vecnames.Push("state");
    scalnames.Push("Error");
    
    an.DefineGraphMesh(2, scalnames, vecnames, "plot.vtk");
    an.PostProcess(0);
    TPZPlasticState<REAL> state = pSD->GetState();
    for (REAL alfa = 1.e-5; alfa<1.e-4; alfa+=1.e-5) {
        state.fAlpha = alfa;
        pSD->SetState(state);
        BuildPlasticSurface(cgrid, pSD);
        an.PostProcess(0);
    }
    delete cgrid;
    delete gmesh;
    
    
}

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD, TPZVec<REAL> &stressresult);

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD)
{
    int ncon = cmesh->NConnects();
    TPZVec<int> computed(ncon,0);
    for (int el=0; el<cmesh->NElements(); el++) {
        TPZCompEl *cel = cmesh->ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        TPZManVector<REAL,3> centerksi(2,0.),xcenter(3,0.);
        gel->CenterPoint(gel->NSides()-1, centerksi);
        gel->X(centerksi, xcenter);
        TPZManVector<REAL> stress(3,0.);
        int matid = ComputeMultiplier(xcenter, pSD,stress);
        cmesh->ElementSolution()(cel->Index(),0) = matid;
        //        gel->SetMaterialId(matid);
        for (int icon=0; icon<gel->NCornerNodes(); icon++)
        {
            TPZConnect &c = cel->Connect(icon);
            TPZGeoNode &gnod = *gel->NodePtr(icon);
            TPZManVector<REAL> co(3,0.),stress(3,0.);
            gnod.GetCoordinates(co);
            
            ComputeMultiplier(co, pSD,stress);
            int seqnum = c.SequenceNumber();
            for (int idf=0; idf<3; idf++) {
                cmesh->Block()(seqnum,0,idf,0) = -co[idf]+stress[idf];
            }
        }
    }
}

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP1> *pSD,TPZVec<REAL> &stressresult)
{
    REAL mult = 1.;
    REAL incr = 1.;
    TPZTensor<REAL> stresstensor,epsilon,stresscenter, epscenter;
    stresstensor.XX() = stress[0];
    stresstensor.YY() = stress[1];
    stresstensor.ZZ() = stress[2];
    stresscenter.XX() = 0.1;
    stresscenter.YY() = 0.1;
    stresscenter.ZZ() = 0.1;
    SANDLERDIMAGGIOSTEP1 *pSDP = dynamic_cast<SANDLERDIMAGGIOSTEP1 *>(pSD);
    pSDP->fER.ComputeDeformation(stresstensor,epsilon);
    pSDP->fER.ComputeDeformation(stresscenter, epscenter);
    TPZTensor<REAL> epsstart(epsilon);
    TPZManVector<REAL,3> phi(2,0.);
    pSD->Phi(epscenter, phi);
    pSD->Phi(epsilon, phi);
    while((phi[0]) > 0. || (phi[1]) > 0.)
    {
        mult *= 0.5;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
    }
    while((phi[0]) < 0. && (phi[1]) < 0.)
    {
        mult *= 2.;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
    }
    REAL tol = mult * 1.e-4;
    mult *= 0.5;
    incr = mult/2.;
    while (incr > tol) {
        mult += incr;
        epsilon = epsstart;
        epsilon.Scale(mult);
        epsilon.Add(epscenter, 1.);
        pSD->Phi(epsilon, phi);
        if ((phi[0]) > 0. || (phi[1]) > 0) {
            mult -= incr;
            incr /= 2.;
        }
    }
    epsilon = epsstart;
    epsilon.Scale(mult);
    epsilon.Add(epscenter, 1.);
    pSD->Phi(epsilon, phi);
    int result = 3;
    if (fabs(phi[0]) < fabs(phi[1])) {
        result = 1;
    }
    else {
        result = 2;
    }
    pSD->fER.Compute(epsilon, stresstensor);
    stressresult[0] = stresstensor.XX();
    stressresult[1] = stresstensor.YY();
    stressresult[2] = stresstensor.ZZ();
    return result;
}
