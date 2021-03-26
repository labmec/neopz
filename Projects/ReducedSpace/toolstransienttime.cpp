//
//  toolstransienttime.cpp
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include <iostream>

#include "toolstransienttime.h"


#include "pzmat1dlin.h"
#include "pzbuildmultiphysicsmesh.h"
#include "pzintel.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "TPZSpStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZJIntegral2D.h"
#include "pzreducedspace.h"
#include "pzbndcond.h"
#include "pzl2projection.h"
#include "tpzmathtools.cpp"
#include "TPZVTKGeoMesh.h"

//Plasticidade
#include "pzelastoplasticanalysis.h"
#include "TPZSandlerDimaggio.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "TPZTensor.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"
#include "TPZMatElastoPlastic2D.h"
//Plasticidade

#include "pzlog.h"
#ifdef PZ_LOG
static TPZLogger logdata("pz.toolstransienttime");
#endif


ToolsTransient::ToolsTransient(){
    fMustStop = true;
    fCouplingMaterial1 = NULL;
    fCouplingMaterial2 = NULL;
    fgmesh = NULL;
    fmeshvec.Resize(0);
    fmphysics = NULL;
    DebugStop();//Nao deveria utilizar este construtor
}

ToolsTransient::ToolsTransient(int pOrder)
{
    fpOrder = pOrder;
    fMustStop = false;
    
    int dim = 2;
    fCouplingMaterial1 = new TPZNLFluidStructure2d(globMultiFisicMatId1,dim,
                                                   globFractInputData.E1(), globFractInputData.Poisson1(), globFractInputData.Visc());
    
    fCouplingMaterial2 = new TPZNLFluidStructure2d(globMultiFisicMatId2,dim,
                                                   globFractInputData.E2(), globFractInputData.Poisson2(), globFractInputData.Visc());
    
    int planestrain = 0;
    fCouplingMaterial1->SetfPlaneProblem(planestrain);
    fCouplingMaterial2->SetfPlaneProblem(planestrain);
    
    fgmesh = NULL;
    fmeshvec.Resize(2);
    fmphysics = NULL;
}

ToolsTransient::~ToolsTransient(){
    
}


void ToolsTransient::RunPlasticity()
{
    //std::map<int,REAL> leakoffMap;
    
    this->Mesh2D();
    
    TPZCompMesh * cmesh = CMeshElastoPlastic(fgmesh, globFractInputData.SigN());
    TPZElastoPlasticAnalysis an(cmesh,std::cout);
    
    this->SolveInitialElastoPlasticity(an, cmesh);
   /* 
    std::cout << "<<<<<<<<<<<<<<<<<<<<<<< 0\n";
    std::stringstream notUsedHere("none.txt");
    REAL KI = ComputeKIPlaneStrain();
    std::cout << "KI = " << KI << " >>>>>>>>>>>>>>>>>>>\n\n\n";
    */
    std::string vtkFile = "pocoplastico.vtk";
    TPZPostProcAnalysis ppanalysis(cmesh);
    ppanalysis.SetStep(0);
    TPZFStructMatrix structmatrix(ppanalysis.Mesh());
    structmatrix.SetNumThreads(8);
    ppanalysis.SetStructuralMatrix(structmatrix);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZStack<std::string> PostProcVars, scalNames, vecNames;
    scalNames.Push("Alpha");
    scalNames.Push("PlasticSqJ2");
    scalNames.Push("PlasticSqJ2El");
    scalNames.Push("POrder");
    
    scalNames.Push("I1Stress");
    scalNames.Push("J2Stress");
    
    vecNames.Push("Displacement");
    vecNames.Push("YieldSurface");
    vecNames.Push("NormalStress");
    vecNames.Push("ShearStress");
    vecNames.Push("NormalStrain");
    vecNames.Push("ShearStrain");
    vecNames.Push("DisplacementMem");
    for (int i=0; i<scalNames.size(); i++)
    {
        PostProcVars.Push(scalNames[i]);
    }
    for (int i=0; i<vecNames.size(); i++)
    {
        PostProcVars.Push(vecNames[i]);
    }
    //
    ppanalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
    //
    ppanalysis.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
    //	
    ppanalysis.TransferSolution();
    ppanalysis.PostProcess(0);// pOrder    
}


void ToolsTransient::Run()
{
    TPZCompMesh * lastMPhysicsCMesh = NULL;
    TPZCompMesh * lastElastReferredCMesh = NULL;
    
    int anCount = 0;
    this->InitializeUncoupledMeshesAttributes();
    this->CMeshMultiphysics();
    TPZAnalysis * an = new TPZAnalysis(fmphysics);
    globFractOutputData.PlotElasticVTK(an, anCount);
    PostprocessPressure();
    PostProcessAcumVolW();
    PostProcessVolLeakoff();
    REAL KI = ComputeKIPlaneStrain();
    globFractOutputData.InsertTKI(globFractInputData.actTime(), KI);//its for output data to txt (Mathematica format)
    
    bool initialElasticKickIsNeeded = true;
    while(fMustStop == false)
    {
        bool propagate = this->SolveSistTransient(an, initialElasticKickIsNeeded);
        initialElasticKickIsNeeded = false;
        anCount = an->GetStep();
        
        if(propagate == true)///Setting new Lfrac and tranferring solution
        {
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
            lastMPhysicsCMesh = fmphysics;
            lastElastReferredCMesh = fmeshvec[0];
        
            REAL newLfrac = globFractInputData.Lf() + globFractInputData.Lmax_edge();
            
            globFractInputData.SetLf(newLfrac);
        
            this->InitializeUncoupledMeshesAttributes();
            this->CMeshMultiphysics();
            this->TransferSolutions(lastMPhysicsCMesh, lastElastReferredCMesh);
            an = new TPZAnalysis(fmphysics);
            
            globFractOutputData.PlotElasticVTK(an, anCount);
            PostProcessAcumVolW();
            PostProcessVolLeakoff();
        }
    }
    
    globFractOutputData.SetQinj1WingAndLfracmax(globFractInputData.Qinj(), globFractInputData.Lf());
    std::ofstream outMath("OutputMathematica.txt");
    globFractOutputData.PrintMathematica(outMath);
}


//------------------------------------------------------------------------------------


void ToolsTransient::InitializeUncoupledMeshesAttributes()
{
    TPZCompMesh * elastReference = ElastCMeshReferenceProcessed();
    
    fmeshvec[0] = this->CMeshReduced(elastReference);
    fmeshvec[1] = this->CMeshPressure();
    fgmesh->ResetReference();
}

TPZCompMesh * ToolsTransient::ElastCMeshReferenceProcessed()
{
    //Principal Geometric Mesh (Lf initial)
    this->Mesh2D();
    
    TPZCompMesh * cmesh_elast = this->CMeshElastic();
    TPZFMatrix<STATE> solutions(cmesh_elast->Solution().Rows(),2);
    
    TPZAnalysis * an = new TPZAnalysis;
    
    an->SetCompMesh(cmesh_elast, true);
    this->SolveInitialElasticity(*an, cmesh_elast);
    
    for(int r = 0; r < cmesh_elast->Solution().Rows(); r++)
    {
        solutions(r,0) = cmesh_elast->Solution()(r,0);
    }
    
    {
        /** Resolvendo um problema modelo de elastica linear para utilizar a
            solucao como espaco de aproximacao do problema nao linear (acoplado) */
        SetSigmaNStripeNum(cmesh_elast,0);
        
        an->SetCompMesh(cmesh_elast, false);
        this->SolveInitialElasticity(*an, cmesh_elast);

        for(int r = 0; r < cmesh_elast->Solution().Rows(); r++)
        {
            solutions(r,1) = cmesh_elast->Solution()(r,0) - solutions(r,0);
        }
    }
    cmesh_elast->LoadSolution(solutions);
    
    return cmesh_elast;
}

void ToolsTransient::Mesh2D()
{
    fgmesh = new TPZGeoMesh;
    
    int ndivV = int(globFractInputData.Lx()/globFractInputData.Lmax_edge() + 0.5);
    int ndivH = int(globFractInputData.Ly()/globFractInputData.Lmax_edge() + 0.5);
    
    int64_t ncols = ndivV + 1;
    int64_t nrows = ndivH + 1;
    int64_t nnodes = nrows*ncols;
    
    fgmesh->NodeVec().Resize(nnodes);
    
    REAL deltadivV = globFractInputData.Lx()/ndivV;
    REAL deltandivH = globFractInputData.Ly()/ndivH;
    
    int64_t nid = 0;
    REAL cracktipDist = globFractInputData.Lf();
    int colCracktip = -1;
    for(int64_t r = 0; r < nrows; r++)
    {
        for(int64_t c = 0; c < ncols; c++)
        {
            REAL x = c*deltadivV;
            REAL y = r*deltandivH;
            REAL dist = fabs(globFractInputData.Lf()-x);
            if(r == 0 && dist < cracktipDist)
            {
                cracktipDist = dist;
                colCracktip = c;
            }
            
            TPZVec<REAL> coord(3,0.);
            coord[0] = x;
            coord[1] = y;
            fgmesh->NodeVec()[r*ncols + c].SetCoord(coord);
            fgmesh->NodeVec()[r*ncols + c].SetNodeId(nid);
            nid++;
        }
    }
    if(colCracktip == 0)
    {
        colCracktip = 1;//fratura minima corresponde aa distancia entre coluna 0 e coluna 1
    }
    
    TPZGeoEl * gel = NULL;
    TPZVec<int64_t> topol(4);
    int64_t indx = 0;
    for(int64_t r = 0; r < nrows-1; r++)
    {
        for(int64_t c = 0; c < ncols-1; c++)
        {
            topol[0] = r*(ncols) + c;
            topol[1] = r*(ncols) + c + 1;
            topol[2] = r*(ncols) + c + 1 + ncols;
            topol[3] = r*(ncols) + c + ncols;
            
            gel = fgmesh->CreateGeoElement(EQuadrilateral, topol, globReservMatId1, indx);
            gel->SetId(indx);
            REAL x = c*deltadivV;
            if((x + 1.E-3) > globFractInputData.Xinterface())
            {
                gel->SetMaterialId(globReservMatId2);
            }
            indx++;
        }
    }
    
    fgmesh->BuildConnectivity();
    
    REAL stripeWidth = globFractInputData.Lf() / globFractInputData.NStripes();
    int64_t nelem = fgmesh->NElements();
    int bcId = globPressureMatId;
    for(int64_t el = 0; el < nelem; el++)
    {
        TPZGeoEl * gel = fgmesh->ElementVec()[el];
        
        //south BC
        TPZGeoElSide sideS(gel,4);
        TPZGeoElSide neighS(sideS.Neighbour());
        if(sideS == neighS)
        {
            if(el < colCracktip)
            {
                TPZGeoEl * bcFrac = gel->CreateBCGeoEl(4,bcId);
                
                //Increasing materialId number with respect with what stripe contains bcFrac
                TPZVec<REAL> centerQsi(bcFrac->Dimension(),0.);
                TPZVec<REAL> centerX(3,0.);
                bcFrac->CenterPoint(bcFrac->NSides()-1, centerQsi);
                bcFrac->X(centerQsi, centerX);
                REAL xCoord = centerX[0];
                int stripeNumber = (int)(xCoord/stripeWidth);

                globFractInputData.InsertBCId_StripeId_ElastId(bcId, stripeNumber, gel->MaterialId());
                
                bcId++;
                ////////
            }
            else
            {
                if(gel->MaterialId() == globReservMatId1)
                {
                    gel->CreateBCGeoEl(4, globDirichletElastMatId1);
                }
                else
                {
                    gel->CreateBCGeoEl(4, globDirichletElastMatId2);
                }
            }
        }
        
        //east BC
        TPZGeoElSide sideE(gel,5);
        TPZGeoElSide neighE(sideE.Neighbour());
        if(sideE == neighE)
        {
            gel->CreateBCGeoEl(5, globDirichletElastMatId2);
        }
        
        //north BC
        TPZGeoElSide sideN(gel,6);
        TPZGeoElSide neighN(sideN.Neighbour());
        if(sideN == neighN)
        {
            if(gel->MaterialId() == globReservMatId1)
            {
                gel->CreateBCGeoEl(6, globDirichletElastMatId1);
            }
            else
            {
                gel->CreateBCGeoEl(6, globDirichletElastMatId2);
            }
        }
        
        //west BC
        TPZGeoElSide sideW(gel,7);
        TPZGeoElSide neighW(sideW.Neighbour());
        if(sideW == neighW)
        {
            gel->CreateBCGeoEl(7, globMixedElastMatId);
        }
    }
    
    topol.Resize(1);
    for(int64_t p = 0; p < ncols; p++)
    {
        topol[0] = p;
        if(p == 0)
        {
            gel = fgmesh->CreateGeoElement(EPoint, topol, globBCfluxIn, indx);
        }
        else if(p == colCracktip)
        {
            gel = fgmesh->CreateGeoElement(EPoint, topol, globCracktip, indx);
        }
        indx++;
    }
    fgmesh->BuildConnectivity();
    
    //#ifdef usingQPoints
    //    TPZGeoElSide pt(gel,0);
    //    TPZGeoElSide ptneigh(pt.Neighbour());
    //    while(pt != ptneigh)
    //    {
    //        if(ptneigh.Element()->HasSubElement() == false)
    //        {
    //            int neighSide = ptneigh.Side();
    //            TPZGeoEl * ptneighEl = TPZChangeEl::ChangeToQuarterPoint(fgmesh, ptneigh.Element()->Id(), neighSide);
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
        nelem = fgmesh->NElements();
        for(int64_t el = 0; el < nelem; el++)
        {
            if(fgmesh->ElementVec()[el]->Dimension() < 1) continue;
            if(fgmesh->ElementVec()[el]->HasSubElement()) continue;
            if(globFractInputData.IsBC(fgmesh->ElementVec()[el]->MaterialId()))
            {
                TPZVec<TPZGeoEl*> sons;
                fgmesh->ElementVec()[el]->Divide(sons);
                continue;
            }
            //else...
            for(int64_t s = 0; s < fgmesh->ElementVec()[el]->NSides(); s++)
            {
                TPZGeoElSide gelside(fgmesh->ElementVec()[el],s);
                TPZGeoElSide neighside(gelside.Neighbour());
                bool refinedAlready = false;
                while(neighside != gelside)
                {
                    if(globFractInputData.IsBC(neighside.Element()->MaterialId()))
                    {
                        TPZVec<TPZGeoEl*> sons;
                        fgmesh->ElementVec()[el]->Divide(sons);
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
    //        nelem = fgmesh->NElements();
    //        for(int el = 0; el < nelem; el++)
    //        {
    //            if(!gmesh->ElementVec()[el]) continue;
    //            if(gmesh->ElementVec()[el]->Dimension() < 1) continue;
    //            if(gmesh->ElementVec()[el]->HasSubElement()) continue;
    //            TPZRefPatternTools::RefineDirectional(gmesh->ElementVec()[el], matDir);
    //        }
    //    }
    //#endif
}

TPZCompMesh * ToolsTransient::CMeshElastic()
{
    /// criar materiais
	int dim = 2;
	
    TPZVec<REAL> force(dim,0.);

    //int planestress = 1;
    int planestrain = 0;
    
    TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                  globFractInputData.E1(),
                                                                  globFractInputData.Poisson1(),
                                                                  globFractInputData.Fx(),
                                                                  globFractInputData.Fy(),
                                                                  planestrain);

    TPZElasticityMaterial * material2 = new TPZElasticityMaterial(globReservMatId2,
                                                                  globFractInputData.E2(),
                                                                  globFractInputData.Poisson2(),
                                                                  globFractInputData.Fx(),
                                                                  globFractInputData.Fy(),
                                                                  planestrain);
    TPZMaterial * mat1(material1);
    TPZMaterial * mat2(material2);
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
    cmesh->SetDefaultOrder(fpOrder);
	cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(mat1);
    cmesh->InsertMaterialObject(mat2);
    
    ///Inserir condicao de contorno
    REAL big = material1->gBigNumber;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    
    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int bcId = it->first;
        int elastId = it->second.second;
        if(elastId == globReservMatId1)
        {//estou no globReservMatId1
            TPZMaterial * BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2);
            cmesh->InsertMaterialObject(BCond11);
        }
        else
        {//estou no globReservMatId2
            TPZMaterial * BCond12 = material2->CreateBC(mat2, bcId, typeNeumann, val1, val2);
            cmesh->InsertMaterialObject(BCond12);
        }
    }
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
    TPZMaterial * BCond22 = material2->CreateBC(mat2, globDirichletElastMatId2, typeDirichlet, val1, val2);
    
    val1(0,0) = big;
    TPZMaterial * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
    
    cmesh->SetAllCreateFunctionsContinuous();
	cmesh->InsertMaterialObject(BCond21);
    cmesh->InsertMaterialObject(BCond22);
	cmesh->InsertMaterialObject(BCond31);
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

void ToolsTransient::SetSigmaNStripeNum(TPZCompMesh * cmeshref, int actStripe)
{
    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int matId = it->first;
        int stripe = it->second.first;
        TPZMaterial * mat = cmeshref->MaterialVec().find(matId)->second;
        if(!mat)
        {
            DebugStop();
        }
        TPZBndCond * bcmat = dynamic_cast<TPZBndCond *>(mat);
        if(!bcmat)
        {
            DebugStop();
        }
        if(stripe < actStripe)
        {
            bcmat->Val2()(1,0) = 0.;
        }
        else if(stripe == actStripe)
        {
            bcmat->Val2()(1,0) = globFractInputData.SigN();
        }
        else if(stripe > actStripe)
        {
            return;
        }
    }
}

TPZCompMeshReferred * ToolsTransient::CMeshReduced(TPZCompMesh *cmeshref){
    /// criar materiais
	int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    int planestrain = 0;
    
    TPZElasticityMaterial * material1 = new TPZElasticityMaterial(globReservMatId1,
                                                                  globFractInputData.E1(),
                                                                  globFractInputData.Poisson1(),
                                                                  globFractInputData.Fx(),
                                                                  globFractInputData.Fy(),
                                                                  planestrain);

    TPZElasticityMaterial * material2 = new TPZElasticityMaterial(globReservMatId2,
                                                                  globFractInputData.E2(),
                                                                  globFractInputData.Poisson2(),
                                                                  globFractInputData.Fx(),
                                                                  globFractInputData.Fy(),
                                                                  planestrain);
    
    material1->SetPreStress(globFractInputData.PreStressXX(), globFractInputData.PreStressYY(), globFractInputData.PreStressXY(), 0.);
    material2->SetPreStress(globFractInputData.PreStressXX(), globFractInputData.PreStressYY(), globFractInputData.PreStressXY(), 0.);
    
    TPZCompMeshReferred * cmeshreferred = new TPZCompMeshReferred(fgmesh);
    
    cmeshreferred->SetDimModel(dim);
    TPZMaterial * mat1(material1);
    TPZMaterial * mat2(material2);
    cmeshreferred->InsertMaterialObject(mat1);
    cmeshreferred->InsertMaterialObject(mat2);
    
    ///Inserir condicao de contorno
    REAL big = material1->gBigNumber;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);

    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int bcId = it->first;
        int elastId = it->second.second;
        if(elastId == globReservMatId1)
        {//estou no globReservMatId1
            TPZMaterial * BCond11 = material1->CreateBC(mat1, bcId, typeNeumann, val1, val2);
            cmeshreferred->InsertMaterialObject(BCond11);
        }
        else
        {//estou no globReservMatId2
            TPZMaterial * BCond12 = material2->CreateBC(mat2, bcId, typeNeumann, val1, val2);
            cmeshreferred->InsertMaterialObject(BCond12);
        }
    }
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond21 = material1->CreateBC(mat1, globDirichletElastMatId1, typeDirichlet, val1, val2);
    TPZMaterial * BCond22 = material2->CreateBC(mat2, globDirichletElastMatId2, typeDirichlet, val1, val2);
    
    val1(0,0) = big;
    TPZMaterial * BCond31 = material1->CreateBC(mat1, globMixedElastMatId, typeMixed, val1, val2);
    
    int numsol = cmeshref->Solution().Cols();
    cmeshreferred->AllocateNewConnect(numsol, 1, 1);
    
	TPZReducedSpace::SetAllCreateFunctionsReducedSpace(cmeshreferred);
    
    cmeshreferred->InsertMaterialObject(BCond21);
    cmeshreferred->InsertMaterialObject(BCond22);
    cmeshreferred->InsertMaterialObject(BCond31);
    
	cmeshreferred->SetDefaultOrder(fpOrder);
    cmeshreferred->SetDimModel(dim);
	
    fgmesh->ResetReference();
	cmeshreferred->AutoBuild();
    cmeshref->AdjustBoundaryElements();
	cmeshref->CleanUpUnconnectedNodes();
    cmeshreferred->LoadReferred(cmeshref);
    
    return cmeshreferred;
}

TPZCompMesh * ToolsTransient::CMeshPressure(){
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(fgmesh);
    cmesh->SetDefaultOrder(fpOrder);
    int dim = 1;
	cmesh->SetDimModel(dim);
    
    /// criar materiais
    TPZFMatrix<REAL> xk(1,1,1.);
    TPZFMatrix<REAL> xc(1,1,0.);
    TPZFMatrix<REAL> xb(1,1,0.);
    TPZFMatrix<REAL> xf(1,1,-2.);
    
    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int bcId = it->first;
        TPZMat1dLin * material = new TPZMat1dLin(bcId);
        material->SetMaterial(xk,xc,xb,xf);
        
        TPZMaterial * mat(material);
    
        cmesh->InsertMaterialObject(mat);
        
        if(it == globFractInputData.GetPressureMatIds_StripeId_ElastId().begin())
        {
            ///Inserir condicao de contorno
            TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
            val2(0,0) = globFractInputData.Qinj();
            TPZMaterial * BCond1 = material->CreateBC(mat, globBCfluxIn, typeNeumann, val1, val2);
            cmesh->InsertMaterialObject(BCond1);
        }

        cmesh->InsertMaterialObject(mat);
    }
    cmesh->SetAllCreateFunctionsContinuous();
	
	//Ajuste da estrutura de dados computacional
	cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
	cmesh->CleanUpUnconnectedNodes();
    
	return cmesh;
}

void ToolsTransient::SolveInitialElasticity(TPZAnalysis &an, TPZCompMesh *Cmesh)
{
	TPZSkylineStructMatrix full(Cmesh); //caso simetrico
	an.SetStructuralMatrix(full);
	TPZStepSolver<REAL> step;
	step.SetDirect(ELDLt); //caso simetrico
	an.SetSolver(step);
	an.Run();
}

void ToolsTransient::CMeshMultiphysics()
{    
    //Creating computational mesh for multiphysic elements
	fgmesh->ResetReference();
	fmphysics = new TPZCompMesh(fgmesh);
    fmphysics->SetDefaultOrder(fpOrder);
    
    TPZMaterial *mat1(fCouplingMaterial1);
    fmphysics->InsertMaterialObject(mat1);
    
    TPZMaterial *mat2(fCouplingMaterial2);
    fmphysics->InsertMaterialObject(mat2);
    
    ///Inserir condicao de contorno
    REAL big = fCouplingMaterial1->gBigNumber;
    
    TPZFMatrix<REAL> val1(3,2,0.), val2(3,1,0.);

    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int bcId = it->first;
        int elastId = it->second.second;
        if(elastId == globReservMatId1)
        {//estou no globReservMatId1
            TPZMaterial * BCond11 = mat1->CreateBC(fCouplingMaterial1, bcId, typeNeumann, val1, val2);
            fmphysics->InsertMaterialObject(BCond11);
        }
        else
        {//estou no globReservMatId2
            TPZMaterial * BCond12 = mat2->CreateBC(fCouplingMaterial2, bcId, typeNeumann, val1, val2);
            fmphysics->InsertMaterialObject(BCond12);
        }
    }
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    TPZMaterial * BCond21 = fCouplingMaterial1->CreateBC(mat1, globDirichletElastMatId1, typeDir_elast, val1, val2);
    TPZMaterial * BCond22 = fCouplingMaterial2->CreateBC(mat2, globDirichletElastMatId2, typeDir_elast, val1, val2);
    
    val1(0,0) = big;
    TPZMaterial * BCond31 = fCouplingMaterial1->CreateBC(mat1, globMixedElastMatId, typeMix_elast, val1, val2);
    
    val2.Redim(3,1);
    val1.Redim(3,2);
    val2(0,0) = globFractInputData.Qinj();
    TPZMaterial * BCond41 = fCouplingMaterial1->CreateBC(mat1, globBCfluxIn, typeNeum_pressure, val1, val2);
    
    fmphysics->SetAllCreateFunctionsMultiphysicElem();
    fmphysics->InsertMaterialObject(BCond21);
    fmphysics->InsertMaterialObject(BCond22);
    fmphysics->InsertMaterialObject(BCond31);
    fmphysics->InsertMaterialObject(BCond41);
    
    fmphysics->AutoBuild();
	fmphysics->AdjustBoundaryElements();
	fmphysics->CleanUpUnconnectedNodes();
    
    // Creating multiphysic elements into mphysics computational mesh
	TPZBuildMultiphysicsMesh::AddElements(fmeshvec, fmphysics);
	TPZBuildMultiphysicsMesh::AddConnects(fmeshvec, fmphysics);
	TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
}


//------------------------------------------------------------------------------------


void ToolsTransient::StiffMatrixLoadVec(TPZAnalysis *an, TPZAutoPointer< TPZMatrix<REAL> > & matK1, TPZFMatrix<REAL> &fvec)
{
	fCouplingMaterial1->SetCurrentState();
    fCouplingMaterial2->SetCurrentState();
    
    TPZFStructMatrix matsk(fmphysics);

	an->SetStructuralMatrix(matsk);
	TPZStepSolver<REAL> step;

	step.SetDirect(ELU);
	an->SetSolver(step);
    
    an->Assemble();
	
    matK1 = an->Solver().Matrix();

	fvec = an->Rhs();
}

void ToolsTransient::TransferSolutions(TPZCompMesh * lastMPhysicsCMesh, TPZCompMesh * lastElastReferredCMesh)
{
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
    TransferElasticSolution(lastElastReferredCMesh);
    TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);

    TransferLeakoff(lastMPhysicsCMesh);
}

void ToolsTransient::TransferElasticSolution(TPZCompMesh * cmeshFrom)
{
#ifdef PZDEBUG
    if(!cmeshFrom)
    {
        DebugStop();
    }
#endif
    
    std::cout << "\n******************** TRANSFERINDO ********************\n\n";
    
    REAL AreaFrom = IntegrateSolution(cmeshFrom, 0);
    
    TPZAutoPointer< TPZFunction<STATE> > func = new TElastSolFunction<STATE>(cmeshFrom);

    //Setting old solution as forcing function on new cmesh (fmeshvec[0])
    TPZMaterial * elastMat = NULL;
    std::map<int,TPZMaterial*>::iterator it = fmeshvec[0]->MaterialVec().find(globReservMatId1);
    if(it != fmeshvec[0]->MaterialVec().end())
    {
        elastMat = it->second;
        elastMat->SetForcingFunction(func);
    }
    else
    {
        it = fmeshvec[0]->MaterialVec().find(globReservMatId2);
        if(it != fmeshvec[0]->MaterialVec().end())
        {
            elastMat = it->second;
            elastMat->SetForcingFunction(func);
        }
        else
        {
            DebugStop();//cade o mardito material???
        }
    }
    
    //////Solving
    TPZAnalysis anTo(fmeshvec[0]);
    TPZSkylineStructMatrix full(fmeshvec[0]);
    anTo.SetStructuralMatrix(full);
    TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);
    anTo.SetSolver(step);
    anTo.Run();    
    anTo.LoadSolution();

    //Restoring original state
    elastMat->SetForcingFunction(NULL);
    
    ///Integral correction
    REAL AreaTo = IntegrateSolution(fmeshvec[0], 0);
    
    if(fabs(AreaTo) > 1.E-18)
    {
        REAL alpha = AreaFrom/AreaTo;
        
        TPZFMatrix<REAL> solutionTo = fmeshvec[0]->Solution();
        for(int r = 0; r < solutionTo.Rows(); r++)
        {
            for(int c = 0; c < solutionTo.Cols(); c++)
            {
                solutionTo(r,c) *= alpha;
            }
        }
        
        fmeshvec[0]->LoadSolution(solutionTo);
    }
}

REAL ToolsTransient::IntegrateSolution(TPZCompMesh * cmesh, int variable)
{
    REAL integral = 0.;
    for(int c = 0; c < cmesh->NElements(); c++)
    {
        TPZCompEl * cel = cmesh->ElementVec()[c];
        if(!cel || globFractInputData.IsBC(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        TPZInterpolationSpace * intel = dynamic_cast<TPZInterpolationSpace*>(cel);
        if(!intel)
        {
            DebugStop();
        }
        TPZVec<REAL> value;
        intel->Integrate(variable, value);
        
        if(value.NElements() == 1)//integrando Leakoff (vlAcum)
        {
            integral += value[0];
        }
        else if(value.NElements() == 2)//Integrando uy (w)
        {
            integral += value[1];
        }
    }
    
    return integral;
}


void ToolsTransient::TransferLeakoff(TPZCompMesh * oldMphysicsCMesh)
{
//    {
//        std::ofstream outAntes("LeakoffANTES.txt");
//        std::map<int,REAL>::iterator it;
//        for(it = globFractInputData.GetLeakoffmap().begin(); it != globFractInputData.GetLeakoffmap().end(); it++)
//        {
//            outAntes << it->second << "\n";
//        }
//    }
    
    TPZCompMesh * cmeshTemp = new TPZCompMesh(fmeshvec[1]->Reference());
    cmeshTemp->SetAllCreateFunctionsDiscontinuous();
    cmeshTemp->AdjustBoundaryElements();
    
    //////L2Projection material
    int dim = 1;
    int pOrder = 0;
    int nsol = 1;
    TPZVec<REAL> solini(nsol,0.);
    
    TPZAutoPointer< TPZFunction<STATE> > func = new TLeakoffFunction<STATE>(oldMphysicsCMesh);
    std::map< int,std::pair<int,int> >::iterator it;
    for(it = globFractInputData.GetPressureMatIds_StripeId_ElastId().begin();
        it != globFractInputData.GetPressureMatIds_StripeId_ElastId().end();
        it++)
    {
        int bcId = it->first;
        TPZL2Projection * materialL2 = new TPZL2Projection(bcId, dim, nsol, solini, pOrder);
        materialL2->SetForcingFunction(func);
    
        //////Inserindo na malha 1D
        cmeshTemp->InsertMaterialObject(materialL2);
    }
    
	cmeshTemp->AutoBuild();
    
    ///////Solving
    TPZAnalysis anTemp(cmeshTemp);
    TPZSkylineStructMatrix full(cmeshTemp);
    anTemp.SetStructuralMatrix(full);
    TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);
    anTemp.SetSolver(step);
    anTemp.Run();
    
    anTemp.LoadSolution();
    
    std::map<int,REAL> newLeakoff;
    for(int cel = 0; cel < cmeshTemp->NElements(); cel++)
    {
        TPZCompEl * compEl = cmeshTemp->ElementVec()[cel];
        TPZGeoEl * geoEl = compEl->Reference();
        
        int gelId = geoEl->Id();
        TPZVec<REAL> center(geoEl->Dimension());
        geoEl->CenterPoint(geoEl->NSides()-1, center);
        
        TPZInterpolationSpace * intpEl = dynamic_cast<TPZInterpolationSpace *>(compEl);
        TPZMaterialData data;
        intpEl->InitMaterialData(data);
        intpEl->ComputeShape(center, data);
        intpEl->ComputeSolution(center, data);
        TPZL2Projection * L2proj = dynamic_cast<TPZL2Projection *>(compEl->Material());
        TPZVec<REAL> Solout(1);
        int var = 1;//TPZL2Projection::ESolution
        L2proj->Solution(data, var, Solout);
        
        REAL vl = Solout[0];
        if(vl < 0.)
        {
            vl = 0.;
        }
        if(geoEl->Neighbour(2).Element()->MaterialId() == globReservMatId1 ||
           geoEl->Neighbour(2).Element()->MaterialId() == globReservMatId2)
        {
            newLeakoff[gelId] = vl;
        }
        else
        {
            DebugStop();
        }
    }
    globFractInputData.GetLeakoffmap() = newLeakoff;
    
//    {
//        std::ofstream outDepois("LeakoffDEPOIS.txt");
//        std::map<int,REAL>::iterator it;
//        for(it = globFractInputData.GetLeakoffmap().begin(); it != globFractInputData.GetLeakoffmap().end(); it++)
//        {
//            outDepois << it->second << "\n";
//        }
//    }
}

void ToolsTransient::MassMatrix(TPZFMatrix<REAL> & Un)
{
    fCouplingMaterial1->SetLastState();
    fCouplingMaterial2->SetLastState();
	TPZSpStructMatrix matsp(fmphysics);
	TPZAutoPointer<TPZGuiInterface> guiInterface;
    matsp.CreateAssemble(Un,guiInterface);
}

bool ToolsTransient::SolveSistTransient(TPZAnalysis *an, bool initialElasticKickIsNeeded)
{
    //CheckConv();
    
	int nrows = an->Solution().Rows();
	TPZFMatrix<REAL> res_total(nrows,1,0.);
    
    TPZFMatrix<REAL> SolIterK = fmphysics->Solution();
	TPZAutoPointer< TPZMatrix<REAL> > matK;
	TPZFMatrix<REAL> fres(fmphysics->NEquations(),1);
    TPZFMatrix<REAL> fmat(fmphysics->NEquations(),1);
    fres.Zero();
    fmat.Zero();
    
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
    
    MassMatrix(fmat);
    
    if(initialElasticKickIsNeeded)
    {
        TPZFMatrix<REAL> chutenewton(fmeshvec[0]->Solution().Rows(), fmeshvec[0]->Solution().Cols(), 1.);
        fmeshvec[0]->LoadSolution(chutenewton);
        TPZBuildMultiphysicsMesh::TransferFromMeshes(fmeshvec, fmphysics);
    }
    
    bool propagate = false;
	while( fMustStop == false && propagate == false )
	{
        fres.Zero();
        StiffMatrixLoadVec(an, matK, fres);
        
        res_total = fres + fmat;
        REAL res = Norm(res_total);
        REAL tol = 1.e-8;
        int maxit = 15;
        int nit = 0;
        
        while(res > tol && nit < maxit) //itercao de Newton
        {
            an->Rhs() = res_total;
            an->Solve();
            an->LoadSolution(SolIterK + an->Solution());
            
            TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
            
            SolIterK = an->Solution();
            
            fres.Zero();
            StiffMatrixLoadVec(an, matK, fres);
            res_total = fres + fmat;

            res = Norm(res_total);
            std::cout << "||res|| = " << res << std::endl;
            nit++;
        }
        
        if(res >= tol)
        {
            std::cout << "\nAtingido o numero maximo de iteracoes, nao convergindo portanto!!!\n";
            std::cout << "||Res|| = " << res << std::endl;
            DebugStop();
        }
        
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
        globFractInputData.UpdateLeakoff(fmeshvec[1]);
        globFractInputData.UpdateActTime();

        PostprocessPressure();
        
        REAL KI = ComputeKIPlaneStrain();
        if(KI > globFractInputData.KIc())
        {//propagou!!!
            globFractInputData.SetMinDeltaT();
            propagate = true;
        }
        else
        {//nao propagou!!!
            globFractInputData.SetNextDeltaT();
            fmat.Zero();
            MassMatrix(fmat);
            
            globFractOutputData.PlotElasticVTK(an);
            PostProcessAcumVolW();
            PostProcessVolLeakoff();
        }
        globFractOutputData.InsertTKI(globFractInputData.actTime(), KI);//its for output data to txt (Mathematica format)
        REAL peteleco = 1.E-8;
        if( globFractInputData.actTime() > (globFractInputData.Ttot() - peteleco) )
        {
            fMustStop = true;
        }
    }
    
    return propagate;
}


REAL ToolsTransient::ComputeKIPlaneStrain()
{
    REAL radius = globFractInputData.Jradius();
    
    REAL KI = -1.;
    
    REAL XcrackTip = -1.;
    TPZGeoMesh * gm = fmeshvec[0]->Reference();
    for(int ell = 0; ell < gm->NElements(); ell++)
    {
        if(gm->ElementVec()[ell] && gm->ElementVec()[ell]->MaterialId() == globCracktip)
        {
            int nodeIndex = gm->ElementVec()[ell]->NodeIndex(0);
            XcrackTip = gm->NodeVec()[nodeIndex].Coord(0);
            break;
        }
    }
    if(XcrackTip < 1.E-3)
    {
        DebugStop();
    }
    TPZVec<REAL> Origin(3,0.);
    Origin[0] = XcrackTip;
    TPZVec<REAL> normalDirection(3,0.);
    normalDirection[2] = 1.;
    
    Path2D * Jpath = new Path2D(fmeshvec[0], Origin, normalDirection, radius);
    JIntegral2D integralJ;
    integralJ.SetPath2D(Jpath);
    
    integralJ.IntegratePath2D();
    REAL computedJ = integralJ.Path()->Jintegral();
    
    REAL young = 0.;
    REAL poisson = 0.;
    
    //<<< Isso nao tah legal qdo o arco da integral J pega os 2 materiais
    if((XcrackTip + 1.E-3) < globFractInputData.Xinterface())
    {
        young = globFractInputData.E1();
        poisson = globFractInputData.Poisson1();
    }
    else
    {
        young = globFractInputData.E2();
        poisson = globFractInputData.Poisson2();
    }
    /////////////////////////////////////////////////////////////////////
    
    if(computedJ < 0.)
    {
        //Estado compressivo!!!
        computedJ = 0.;
    }
    KI = sqrt( (computedJ*young) / (1. - poisson*poisson) );
    
    std::cout << "\n>>>>>>>>>>> KI = " << KI << "\n";
    
    return KI;
}

void ToolsTransient::PostprocessPressure()
{
    std::map<REAL,REAL> pos_pressure;

    for(int i = 0;  i < fmeshvec[1]->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel = fmeshvec[1]->ElementVec()[i];
        TPZInterpolatedElement * sp = dynamic_cast <TPZInterpolatedElement*>(cel);
        if(!sp)
        {
            continue;
        }
        TPZMaterialData data;
        sp->InitMaterialData(data);

        for(int qsiPt = -1; qsiPt <= +1; qsiPt++)
        {
            TPZVec<REAL> qsi(1,0.), Xqsi(3,0.);
            qsi[0] = qsiPt;
            cel->Reference()->X(qsi,Xqsi);
            
            sp->ComputeShape(qsi, data);
            sp->ComputeSolution(qsi, data);

            REAL pos = Xqsi[0];
            REAL press = data.sol[0][0];
            pos_pressure[pos] = press;
        }
    }
    
    globFractOutputData.InsertTposP(globFractInputData.actTime(), pos_pressure);
}

void ToolsTransient::PostProcessAcumVolW()
{
    REAL wAcum = 2. * IntegrateSolution(fmeshvec[0], 0);//aqui jah sao consideradas as 2 faces da fratura
    globFractOutputData.InsertTAcumVolW(globFractInputData.actTime(), wAcum);
}

void ToolsTransient::PostProcessVolLeakoff()
{
    //volume por trecho
    TPZAutoPointer< TPZFunction<STATE> > func = new TLeakoffFunction<STATE>(fmphysics);
    
    REAL vlAcum = 0.;

    for(int i = 0;  i < fmeshvec[1]->ElementVec().NElements(); i++)
    {
        TPZCompEl * cel = fmeshvec[1]->ElementVec()[i];
        if(globFractInputData.IsBC(cel->Reference()->MaterialId()) == false)
        {
            continue;
        }
        
#ifdef PZDEBUG
        if(!cel || cel->Reference()->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        cel->Material()->SetForcingFunction(func);
        TPZGeoEl * fluidGel = cel->Reference();
        
#ifdef PZDEBUG
        if(!fluidGel)
        {
            DebugStop();
        }
#endif
        
        TPZVec<REAL> qsi(1,0.), xLeft(3,0.), xMiddle(3,0.), xRight(3,0.), leakoff(1,0);
        
        qsi[0] = -1.;
        fluidGel->X(qsi,xLeft);
        REAL posLeft = ((int)(xLeft[0]*1000.))/1000.;
        
        qsi[0] = +1.;
        fluidGel->X(qsi,xRight);
        REAL posRight = ((int)(xRight[0]*1000.))/1000.;
        
        REAL lengthElement = xRight[0] - xLeft[0];
        
        xMiddle[0] = (xLeft[0] + xRight[0])/2.;
        func->Execute(xMiddle, leakoff);
        
        //soh para nao coincidir pontos no mapa
        xLeft[0] = xMiddle[0] - 0.99*lengthElement/2.;
        xRight[0] = xMiddle[0] + 0.99*lengthElement/2.;

        REAL lengthLeakoff = leakoff[0];
        globFractOutputData.InsertTposVolLeakoff(globFractInputData.actTime(), posRight, lengthLeakoff);
        globFractOutputData.InsertTposVolLeakoff(globFractInputData.actTime(), posLeft, lengthLeakoff);
        
        cel->Material()->SetForcingFunction(NULL);
        
        //volume Acumulado
        vlAcum += 2. * (lengthLeakoff * lengthElement);//sao duas faces da fratura, lembra?
    }
    
    globFractOutputData.InsertTAcumVolLeakoff(globFractInputData.actTime(), vlAcum);
}

void ToolsTransient::CheckConv()
{    
	int64_t neq = fmphysics->NEquations();
    int nsteps = 10;
    
    TPZFMatrix<REAL> xIni(neq,1);
    for(int64_t i = 0; i < xIni.Rows(); i++)
    {
        REAL val = (double)(rand())*(1.e-10);
        xIni(i,0) = val;
    }
    TPZAnalysis *an = new TPZAnalysis(fmphysics);
    an->LoadSolution(xIni);
    TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
    
    TPZFMatrix<REAL> actX = xIni;
    
    TPZAutoPointer< TPZMatrix<REAL> > fL_xIni;
    TPZFMatrix<REAL> f_xIni(neq,1);
    
    StiffMatrixLoadVec(an, fL_xIni, f_xIni);
    if(fL_xIni->Rows() != neq || fL_xIni->Cols() != neq || fL_xIni->IsDecomposed())
    {
        DebugStop();
    }
    
    TPZFMatrix<REAL> fAprox_x(neq,1);
    TPZFMatrix<REAL> fExato_x(neq,1);
    
    TPZFMatrix<REAL> errorVec(neq,1,0.);
	TPZFMatrix<REAL> errorNorm(nsteps,1,0.);
    
    
    TPZAutoPointer< TPZMatrix<REAL> > fLtemp;
    TPZFMatrix<REAL> dFx(neq,1);
    
    TPZVec<REAL> deltaX(neq,0.01), alphas(nsteps);
    double alpha;
    
    std::stringstream exatoSS, aproxSS;
    exatoSS << "exato={";
    aproxSS << "aprox={";
    for(int i = 0; i < nsteps; i++)
    {
        alpha = (i+1)/10.;
        alphas[i] = alpha;
        
        ///Fx aproximado
        dFx.Zero();
        for(int64_t r = 0; r < neq; r++)
        {
            for(int64_t c = 0; c < neq; c++)
            {
                dFx(r,0) +=  (-1.) * fL_xIni->GetVal(r,c) * (alpha * deltaX[c]); // (-1) porque fLini = -D[res,sol]
            }
        }
        fAprox_x = f_xIni + dFx;
        
        int wantToSeeRow = 0;
        REAL aproxSol = fAprox_x(wantToSeeRow,0);
        std::cout << "Aprox : " << aproxSol << std::endl;
        {
            aproxSS << aproxSol;
            if(i < nsteps-1)
            {
                aproxSS << ",";
            }
        }
        
        ///Fx exato
        for(int64_t r = 0; r < neq; r++)
        {
            actX(r,0) = xIni(r,0) + (alpha * deltaX[r]);
        }
        an->LoadSolution(actX);
        TPZBuildMultiphysicsMesh::TransferFromMultiPhysics(fmeshvec, fmphysics);
        
        fExato_x.Zero();
        if(fLtemp) fLtemp->Zero();
        StiffMatrixLoadVec(an, fLtemp, fExato_x);
        
        REAL exatoSol = fExato_x(wantToSeeRow,0);
        std::cout << "Exato : " << exatoSol << std::endl;
        {
            exatoSS << exatoSol;
            if(i < nsteps-1)
            {
                exatoSS << ",";
            }
        }
        
        ///Erro
        errorVec.Zero();
        for(int64_t r = 0; r < neq; r++)
        {
            errorVec(r,0) = fExato_x(r,0) - fAprox_x(r,0);
        }
        
        ///Norma do erro
        double XDiffNorm = Norm(errorVec);
        errorNorm(i,0) = XDiffNorm;
    }
    aproxSS << "};";
    exatoSS << "};";
    std::cout << aproxSS.str() << std::endl;
    std::cout << exatoSS.str() << std::endl;
    std::cout << "Show[ListPlot[aprox, Joined -> True, PlotStyle -> Red],ListPlot[exato, Joined -> True]]\n";
    
    std::cout << "Convergence Order:\n";
    for(int j = 1; j < nsteps; j++)
    {
        std::cout << ( log(errorNorm(j,0)) - log(errorNorm(j-1,0)) )/( log(alphas[j]) - log(alphas[j-1]) ) << "\n";
    }
}


//---------------------------------------------------------------


template<class TVar>
TElastSolFunction<TVar>::TElastSolFunction()
{
    fIniElIndex = 0;
    fcmesh = NULL;
}

template<class TVar>
TElastSolFunction<TVar>::TElastSolFunction(TPZCompMesh * cmesh)
{
    fIniElIndex = 0;
    this->fcmesh = cmesh;
}

template<class TVar>
TElastSolFunction<TVar>::~TElastSolFunction()
{
    fIniElIndex = 0;
    fcmesh = NULL;
}

template<class TVar>
void TElastSolFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f)
{
    TPZVec<REAL> xcp(x);
    TPZVec<REAL> qsi2D(2,0.);
    TPZGeoEl * gel = fcmesh->Reference()->FindElement(xcp, qsi2D, fIniElIndex, 2);

    if(!gel)
    {
        DebugStop();
    }
    
    TPZElasticityMaterial * elast = dynamic_cast<TPZElasticityMaterial *>(gel->Reference()->Material());
    if(!elast)
    {
        DebugStop();
    }
    
    f.Resize(2);
    
    TPZReducedSpace * intpEl = dynamic_cast<TPZReducedSpace *>(gel->Reference());
    if(!intpEl)
    {
        DebugStop();
    }
    TPZMaterialData data;
    intpEl->InitMaterialData(data);
    
    intpEl->ComputeShape(qsi2D, data);
    intpEl->ComputeSolution(qsi2D, data);
    
    int var = 0;
    elast->Solution(data, var, f);
}

template<class TVar>
void TElastSolFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
{
    DebugStop();
}

template<class TVar>
int TElastSolFunction<TVar>::NFunctions() const
{
    return 2;
}

template<class TVar>
int TElastSolFunction<TVar>::PolynomialOrder() const
{
    return fcmesh->GetDefaultOrder();
}


//---------------------------------------------------------------


template<class TVar>
TLeakoffFunction<TVar>::TLeakoffFunction()
{
    DebugStop();//use o outro construtor
}

template<class TVar>
TLeakoffFunction<TVar>::TLeakoffFunction(TPZCompMesh * cmesh)
{
    this->fIniElIndex = 0;
    this->fcmesh = cmesh;
}

template<class TVar>
TLeakoffFunction<TVar>::~TLeakoffFunction()
{
    fIniElIndex = 0;
    fcmesh = NULL;
}

template<class TVar>
void TLeakoffFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f)
{
    TPZVec<REAL> xcp(x);
    TPZVec<REAL> qsi1D(1,0.);
    TPZGeoEl * gel = fcmesh->Reference()->FindElement(xcp, qsi1D, fIniElIndex, 1);
    if(!gel)
    {
        DebugStop();
    }
    
    if(globFractInputData.IsBC(gel->MaterialId()) == false)
    {
        TPZGeoElSide gelSide(gel->NSides()-1);
        TPZGeoElSide neighSide(gelSide.Neighbour());
        while(neighSide != gelSide)
        {
            if(globFractInputData.IsBC(neighSide.Element()->MaterialId()))
            {
                gel = neighSide.Element();
                break;
            }
            neighSide = neighSide.Neighbour();
        }
    }
    if(globFractInputData.IsBC(gel->MaterialId()) == false)
    {
        f[0] = 0.;
        return;
    }
    
    f.Resize(1);
    int elId = gel->Id();
    std::map<int,REAL>::iterator it = globFractInputData.GetLeakoffmap().find(elId);
    if(it == globFractInputData.GetLeakoffmap().end())
    {
        f[0] = 0.;
    }
    else
    {
        f[0] = it->second;
    }
}

template<class TVar>
void TLeakoffFunction<TVar>::Execute(const TPZVec<REAL> &x, TPZVec<TVar> &f, TPZFMatrix<TVar> &df)
{
    DebugStop();
}

template<class TVar>
int TLeakoffFunction<TVar>::NFunctions() const
{
    return 1;
}

template<class TVar>
int TLeakoffFunction<TVar>::PolynomialOrder() const
{
    return 0;
}

void ToolsTransient::SolveInitialElastoPlasticity(TPZElastoPlasticAnalysis &analysis, TPZCompMesh *Cmesh)
{
	TPZSkylineStructMatrix full(Cmesh);
    full.SetNumThreads(0);
	analysis.SetStructuralMatrix(full);
    
	TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);
	analysis.SetSolver(step);
    
    int NumIter = 2;
    bool linesearch = true;
    bool checkconv = false;
    
    analysis.IterativeProcess(cout, 1.e-6, NumIter, linesearch, checkconv);
    
    analysis.AcceptSolution();
}

TPZCompMesh * ToolsTransient::CMeshElastoPlastic(TPZGeoMesh *gmesh, REAL SigmaN)
{
    /// criar materiais
    int dim = 2;
    
    TPZVec<REAL> force(dim,0.);
    
    //int planestress = 1;
    int planestrain = 1;
    
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> SD;
    REAL poisson = 0.203;
    REAL elast = 29269.;
    REAL A = 152.54;
    REAL B = 0.0015489;
    REAL C = 146.29;
    REAL R = 0.91969;
    REAL D = 0.018768;
    REAL W = 0.006605;
    SD.SetUp(poisson, elast, A, B, C, R, D, W);
    SD.SetResidualTolerance(1.e-10);
    SD.fIntegrTol = 10.;
    
    TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > * PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > (globReservMatId1,planestrain);
    
    TPZMaterial * mat(PlasticSD);
    PlasticSD->SetPlasticityModel(SD);
    
    TPZMatWithMem<TPZElastoPlasticMem> * pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (mat);
    
    ///criar malha computacional
    TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
    cmesh->SetDefaultOrder(fpOrder);
    cmesh->SetDimModel(dim);
    cmesh->InsertMaterialObject(pMatWithMem);
    
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(cmesh);
    
    
    ///Inserir condicao de contorno
    REAL big = mat->gBigNumber;
    
    TPZFMatrix<REAL> val1(2,2,0.), val2(2,1,0.);
    val2(1,0) = SigmaN;
    TPZMaterial * BCond1 = PlasticSD->CreateBC(pMatWithMem, globPressureMatId, typeNeumann, val1, val2);///NATHAN, falta incrementar stripenumber!!!
    
    val1.Redim(2,2);
    val2.Redim(2,1);
    TPZMaterial * BCond2 = PlasticSD->CreateBC(pMatWithMem, globDirichletElastMatId1, typeDirichlet, val1, val2);///NATHAN, agora temos 2 faixas de mat elastico!!!
    
    val1(0,0) = big;
    TPZMaterial * BCond3 = PlasticSD->CreateBC(pMatWithMem, globMixedElastMatId, typeMixed, val1, val2);
    
    //cmesh->SetAllCreateFunctionsContinuous();
    cmesh->InsertMaterialObject(BCond1);
    cmesh->InsertMaterialObject(BCond2);
    cmesh->InsertMaterialObject(BCond3);
    
    //Ajuste da estrutura de dados computacional
    cmesh->AutoBuild();
    cmesh->AdjustBoundaryElements();
    cmesh->CleanUpUnconnectedNodes();
    
    return cmesh;
}
