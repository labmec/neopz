//
//  File.cpp
//  PZ
//
//  Created by Agnaldo Farias on 7/31/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//


#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzreferredcompel.h"
#include "pzcompel.h"
#include "pzbndcond.h"
#include "TPZInterfaceEl.h"
#include "pzinterpolationspace.h"

#include "pzpoisson3d.h"
#include "pzelasmat.h"
#include "pzmat1dlin.h"
#include "pzelastpressure.h"
#include "pznlfluidstructure2d.h"

#include "tpzgeoelrefpattern.h"
#include "TPZGeoLinear.h"
#include "tpztriangle.h"
#include "pzgeoquad.h"
#include "pzgeopoint.h"

#include "pzanalysis.h"
#include "pznonlinanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstrmatrix.h"
#include "pzstepsolver.h"
#include "pzreducedspace.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"

#include "pzbuildmultiphysicsmesh.h"

#include "TPZVTKGeoMesh.h"

#include <iostream>
#include <math.h>

#include <fstream>
#include <sstream>

#include "toolstransienttime.h"

#include "tpzchangeel.h"

using namespace std;

int const reservMatId   = 1; //elastic
int const pressureMatId = 2; //pressure

int const dirichletElastMatId = -1;
int const neumannElastMatId   = -2;
int const blockedXElastMatId  = -3;

int const bcfluxIn  = -10; //bc pressure
int const bcfluxOut = -20; //bc pressure

int const dirichlet = 0;
int const neumann   = 1;
int const mixed     = 2;

int const neum_elast    = 10;
int const dir_elast     = 11;
int const mix_elast     = 20;
int const neum_pressure = 21;
int const dir_pressure  = 22;


REAL const Pi = 4.*atan(1.);

TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax);

TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh, int pOrder);
TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder);
TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh, int pOrder);

TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial);
void SaidaPressao(TPZVec<TPZCompMesh *> meshvec, TPZCompMesh* mphysics);

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh);

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.reducedspace.data"));
#endif

//Dimensions:
//const REAL Lx = 2.;
//const REAL Ly = 2.;
//const REAL Lf = 1.;
//const REAL Hf = 1.;
////Elastic properties:
//const REAL ED = 0.3E5;
//const REAL nu = 0.25;
////Fluid property:
//const REAL visc = 0.001;
////BCs:
//const REAL sigN = 0.3E5;/// <<< sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
//const REAL Qinj  = -1.E-3;
////time:
//const REAL Ttot = 1000.;
//const REAL deltaT = 10.;
////Leakoff:
//const REAL Cl = 0.0007;
//const REAL Pe = 1.;
//const REAL Pref = 2000.;
//const REAL vsp = 0.00001;
//



 //Dimensions:
 const REAL Lx = 400.;
 const REAL Ly = 400.;
 REAL Lf = 50.;
 const REAL Hf = 1.;
 //Elastic properties:
 const REAL ED = 3.9E4;//MPa
 const REAL nu = 0.25;
 //Fluid property:
 const REAL visc = 0.001E-6;//MPa.s
 //BCs:
 const REAL sigN = 61.5;/// <<< sigma.n no problema elastico que servira de espaco de aproximacao para o elastico multifisico
 const REAL Qinj  = -50./Hf;///vazao de 1 asa de fratura dividido pela altura da fratura
 //time:
 const REAL Ttot = 50.;//sem bug -> Ttot = 50.
 const REAL nsteps = 20.;//sem bug -> nsteps = 20.
 const REAL deltaT = Ttot/nsteps;
 //Leakoff:
 const REAL Cl = 0.01;
 const REAL Pe = 10.;//MPa
 const REAL SigmaConf = 11.;//MPa
 const REAL Pref = 60000.;//MPa
 const REAL vsp = 0.001;
 //
 


int main(int argc, char *argv[])
{
#ifdef LOG4CXX
    InitializePZLOG();
#endif
    
#ifdef LOG4CXX
//    std::stringstream str;
//    str << "TESTANDO!!!";
//    LOGPZ_DEBUG(logger,str.str());    
#endif
    
    int p = 2;
	//primeira malha
	
	// geometric mesh (initial)
    TPZGeoMesh * gmesh = PlaneMesh(Lf, Lx, Ly, 10.);
    
    //computational mesh elastic
/** >>> Resolvendo um problema modelo de elastica linear para utilizar a solucao 
    como espaco de aproximacao do problema nao linear (acoplado) */
	TPZCompMesh * cmesh_elast = CMeshElastic(gmesh, p);
    TPZAnalysis an1(cmesh_elast);
    MySolve(an1, cmesh_elast);
    
/** >>> Passando a malha computacional jah processada do problema modelo elastico para a malha do problema elastico que serah acoplado.
    Esta malha que foi passada servirah como espaco de aproximacao da malha do problema elastico que serah acoplado */
    TPZCompMeshReferred * cmesh_referred = CMeshReduced(gmesh, cmesh_elast, p);
    cmesh_referred->ComputeNodElCon();
    TPZFStructMatrix fstr(cmesh_referred);
    TPZFMatrix<STATE> rhs;//(1);
    TPZAutoPointer<TPZMatrix<STATE> > strmat = fstr.CreateAssemble(rhs,NULL);
	
/** >>> Criando a malha computacional para o problema de fluxo de fluido */
    TPZCompMesh * cmesh_pressure = CMeshPressure(gmesh, p);
    
/** >>> Problema multifisico : acoplamento das 2 malhas computacionais anteriores (elastica e fluxo) */
    //multiphysic mesh
    TPZVec<TPZCompMesh *> meshvec(2);
	meshvec[0] = cmesh_referred;
	meshvec[1] = cmesh_pressure;
    
    gmesh->ResetReference();
	TPZNLFluidStructure2d * mymaterial = NULL;
/** >>> Serah utilizado o material criado (pznlfluidstructure2D do Agnaldo) */
    TPZCompMesh * mphysics = MalhaCompMultphysics(gmesh, meshvec, mymaterial);
    mphysics->SetDefaultOrder(p);
    
    mymaterial->SetTimeStep(deltaT);
    REAL maxTime = Ttot;
    
/** >>> Metodo de resolucao de problema transient */
    TPZAnalysis *an = new TPZAnalysis(mphysics);
    TPZFMatrix<REAL> InitialSolution = an->Solution();
    ToolsTransient::SolveSistTransient(deltaT, maxTime, InitialSolution, an, mymaterial, meshvec, mphysics);
    
    return 0;
}


TPZGeoMesh * PlaneMesh(REAL lf, REAL ldom, REAL hdom, REAL lmax)
{
    TPZGeoMesh * gmesh = new TPZGeoMesh;
    
    int ndivV = int(ldom/lmax + 0.5);
    int ndivH = int(hdom/lmax + 0.5);
    
    int ncols = ndivV + 1;
    int nrows = ndivH + 1;
    int nnodes = nrows*ncols;
    
    gmesh->NodeVec().Resize(nnodes);
    
    REAL deltadivV = ldom/ndivV;
    REAL deltandivH = hdom/ndivH;
    
    int nid = 0;
    REAL cracktipDist = lf;
    int colCracktip = -1;
    for(int r = 0; r < nrows; r++)
    {
        for(int c = 0; c < ncols; c++)
        {
            REAL x = c*deltadivV;
            REAL y = r*deltandivH;
            REAL dist = fabs(lf-x);
            if(dist < cracktipDist)
            {
                cracktipDist = dist;
                colCracktip = c;
            }
            
            TPZVec<REAL> coord(3,0.);
            coord[0] = x;
            coord[1] = y;
            gmesh->NodeVec()[r*ncols + c].SetCoord(coord);
            gmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
            nid++;
        }
    }
    if(colCracktip == 0)
    {
        colCracktip = 1;//fratura minima corresponde aa distancia entre coluna 0 e coluna 1
    }
    
    TPZGeoEl * gel = NULL;
    TPZVec<int> topol(4);
    int indx = 0;
    for(int r = 0; r < nrows-1; r++)
    {
        for(int c = 0; c < ncols-1; c++)
        {
            topol[0] = r*(ncols) + c;
            topol[1] = r*(ncols) + c + 1;
            topol[2] = r*(ncols) + c + 1 + ncols;
            topol[3] = r*(ncols) + c + ncols;
            
            gel = gmesh->CreateGeoElement(EQuadrilateral, topol, reservMatId, indx);
            gel->SetId(indx);
            indx++;
        }
    }
    
    gmesh->BuildConnectivity();
    
    int nelem = gmesh->NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = gmesh->ElementVec()[el];
        
        //south BC
        TPZGeoElSide sideS(gel,4);
        TPZGeoElSide neighS(sideS.Neighbour());
        if(sideS == neighS)
        {
            if(el < colCracktip)
            {
                gel->CreateBCGeoEl(4, pressureMatId);
            }
            else
            {
                gel->CreateBCGeoEl(4, dirichletElastMatId);
            }
        }
        
        //east BC
        TPZGeoElSide sideE(gel,5);
        TPZGeoElSide neighE(sideE.Neighbour());
        if(sideE == neighE)
        {
            gel->CreateBCGeoEl(5, dirichletElastMatId);
        }
        
        //north BC
        TPZGeoElSide sideN(gel,6);
        TPZGeoElSide neighN(sideN.Neighbour());
        if(sideN == neighN)
        {
            gel->CreateBCGeoEl(6, dirichletElastMatId);
        }
        
        //west BC
        TPZGeoElSide sideW(gel,7);
        TPZGeoElSide neighW(sideW.Neighbour());
        if(sideW == neighW)
        {
            gel->CreateBCGeoEl(7, blockedXElastMatId);
        }
    }
    
    topol.Resize(1);
    topol[0] = 0;
    gel = gmesh->CreateGeoElement(EPoint, topol, bcfluxIn, indx);
    indx++;
    
    topol[0] = colCracktip;
    gel = gmesh->CreateGeoElement(EPoint, topol, bcfluxOut, indx);
    
    gmesh->BuildConnectivity();
    
//#ifdef usingQPoints
//    TPZGeoElSide pt(gel,0);
//    TPZGeoElSide ptneigh(pt.Neighbour());
//    while(pt != ptneigh)
//    {
//        if(ptneigh.Element()->HasSubElement() == false)
//        {
//            int neighSide = ptneigh.Side();
//            TPZGeoEl * ptneighEl = TPZChangeEl::ChangeToQuarterPoint(gmesh, ptneigh.Element()->Id(), neighSide);
//            ptneigh = ptneighEl->Neighbour(neighSide);
//        }
//        else
//        {
//            ptneigh = ptneigh.Neighbour();
//        }
//    }
//#endif
    
    int nrefUnif = 3;
    for(int ref = 0; ref < nrefUnif; ref++)
    {
        nelem = gmesh->NElements();
        for(int el = 0; el < nelem; el++)
        {
            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
            if(gmesh->ElementVec()[el]->MaterialId() == pressureMatId)
            {
                TPZVec<TPZGeoEl*> sons;
                gmesh->ElementVec()[el]->Divide(sons);
                continue;
            }
            //else...
            for(int s = 0; s < gmesh->ElementVec()[el]->NSides(); s++)
            {
                TPZGeoElSide gelside(gmesh->ElementVec()[el],s);
                TPZGeoElSide neighside(gelside.Neighbour());
                bool refinedAlready = false;
                while(neighside != gelside)
                {
                    if(neighside.Element()->MaterialId() == pressureMatId)
                    {
                        TPZVec<TPZGeoEl*> sons;
                        gmesh->ElementVec()[el]->Divide(sons);
                        refinedAlready = true;
                        break;
                    }
                    neighside = neighside.Neighbour();
                }
                if(refinedAlready == true)
                {
                    break;
                }
            }
        }
    }
    
//#ifdef usingRefdir
//    std::set<int> matDir;
//    //matDir.insert(__2DfractureMat_inside);
//    matDir.insert(bcfluxOut);
//    int nrefDir = 1;
//    for(int ref = 0; ref < nrefDir; ref++)
//    {
//        nelem = gmesh->NElements();
//        for(int el = 0; el < nelem; el++)
//        {
//            if(!gmesh->ElementVec()[el]) continue;
//            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
//            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
//            TPZRefPatternTools::RefineDirectional(gmesh->ElementVec()[el], matDir);
//        }
//    }
//#endif
    
    return gmesh;
}

TPZCompMesh * CMeshElastic(TPZGeoMesh *gmesh, int pOrder)
{
    /// criar materiais
	int dim = 2;
	
    TPZVec<REAL> force(dim,0.);
    REAL E = ED;
	REAL poisson = nu;
    //int planestress = 1;
    int planestrain = 0;
    
    TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(reservMatId, E, poisson, force[0], force[1], planestrain);

    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(1,0) = sigN;
    TPZMaterial * BCond1 = material->CreateBC(mat, pressureMatId, neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond2 = material->CreateBC(mat, dirichletElastMatId, dirichlet, val1, val2);
     
    val1(0,0) = big;
    TPZMaterial * BCond3 = material->CreateBC(mat, blockedXElastMatId, mixed, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(BCond1);
	cmesh->InsertMaterialObject(BCond2);
	cmesh->InsertMaterialObject(BCond3);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
        
	return cmesh;
}

TPZCompMeshReferred * CMeshReduced(TPZGeoMesh *gmesh, TPZCompMesh *cmesh, int pOrder){
    
    /// criar materiais
	int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    REAL E = ED;
	REAL poisson = nu;
    //int planestress = 1;
    int planestrain = 0;
    
    TPZElasticityMaterial *material;
	material = new TPZElasticityMaterial(reservMatId, E, poisson, force[0], force[1], planestrain); 
	material->NStateVariables();
    
    
    TPZCompMeshReferred *cmeshreferred = new TPZCompMeshReferred(gmesh);
    
    cmeshreferred->SetDimModel(dim);
    TPZMaterial * mat(material);
    cmeshreferred->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    REAL big = material->gBigNumber;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    TPZMaterial * BCond1 = material->CreateBC(mat, pressureMatId, neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond2 = material->CreateBC(mat, dirichletElastMatId, dirichlet, val1, val2);
    
    val1(0,0) = big;
    TPZMaterial * BCond3 = material->CreateBC(mat, blockedXElastMatId, mixed, val1, val2);

    int numsol = cmesh->Solution().Cols();
    cmeshreferred->AllocateNewConnect(numsol, 1, 1);
    
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
    
    cmeshreferred->InsertMaterialObject(BCond1);
    cmeshreferred->InsertMaterialObject(BCond2);
    cmeshreferred->InsertMaterialObject(BCond3);
    
	cmeshreferred->SetDefaultOrder(pOrder);
    cmeshreferred->SetDimModel(dim);
	
    gmesh->ResetReference();
	cmeshreferred->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    cmeshreferred->LoadReferred(cmesh);
   
    return cmeshreferred;
}

TPZCompMesh * CMeshPressure(TPZGeoMesh *gmesh, int pOrder){
    
    /// criar materiais
	int dim = 1;
	
    TPZMat1dLin *material;
	material = new TPZMat1dLin(pressureMatId);
    
    TPZFMatrix<REAL> xk(1,1,1.);
    TPZFMatrix<REAL> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<REAL> xf(1,1,-2.);
    material->SetMaterial(xk,xc,xb,xf);
    
    TPZMaterial * mat(material);
    material->NStateVariables();
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(pOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat);
    
    ///Inserir condicao de contorno
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(0,0) = Qinj;
    TPZMaterial * BCond1 = material->CreateBC(mat, bcfluxIn, neumann, val1, val2);
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond2 = material->CreateBC(mat, bcfluxOut, dirichlet, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(mat);
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;    
}

void MySolve(TPZAnalysis &an, TPZCompMesh *Cmesh)
{			
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	an.SetSolver(step);
	an.Run();
}

TPZCompMesh * MalhaCompMultphysics(TPZGeoMesh * gmesh, TPZVec<TPZCompMesh *> meshvec, TPZNLFluidStructure2d * &mymaterial){
    
    //Creating computational mesh for multiphysic elements
	gmesh->ResetReference();
	TPZCompMesh *mphysics = new TPZCompMesh(gmesh);
    
    int MatId=1;
    int dim =2;
    mymaterial = new TPZNLFluidStructure2d(MatId,dim);
    
    //data elasticity
    REAL fx = 0.;
    REAL fy = 0.;
    
    //int planestress = 1;
    int planestrain = 0;
    mymaterial->SetElasticParameters(ED, nu, fx, fy);
    mymaterial->SetfPlaneProblem(planestrain);

    mymaterial->SetParameters(Hf, Lf, visc, Qinj, Cl, Pe, SigmaConf, Pref, vsp);
    
    TPZMaterial *mat(mymaterial);
    mphysics->InsertMaterialObject(mat);
    
    mymaterial->NStateVariables();
    ///Inserir condicao de contorno
    REAL big = mymaterial->gBigNumber;
    
    TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);
    val2(1,0) = sigN;
    TPZMaterial * BCond1 = mymaterial->CreateBC(mat, pressureMatId, neumann, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    TPZMaterial * BCond2 = mymaterial->CreateBC(mat, dirichletElastMatId, dir_elast, val1, val2);
    
    val1(0,0) = big;
    TPZMaterial * BCond3 = mymaterial->CreateBC(mat, blockedXElastMatId, mix_elast, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    val2(2,0) = Qinj;
    TPZMaterial * BCond4 = mymaterial->CreateBC(mat, bcfluxIn, neum_pressure, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    TPZMaterial * BCond5 = mymaterial->CreateBC(mat, bcfluxOut, neum_pressure, val1, val2);
    
    mphysics->SetAllCreateFunctionsMultiphysicElem();
    mphysics->InsertMaterialObject(BCond1);
    mphysics->InsertMaterialObject(BCond2);
    mphysics->InsertMaterialObject(BCond3);
    mphysics->InsertMaterialObject(BCond4);
    mphysics->InsertMaterialObject(BCond5);

    mphysics->AutoBuild();
	mphysics->AdjustBoundaryElements();
	mphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(meshvec, mphysics);
	TPZBuildMultiphysicsMesh::AddConnects(meshvec,mphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(meshvec, mphysics);

    return mphysics;
}

