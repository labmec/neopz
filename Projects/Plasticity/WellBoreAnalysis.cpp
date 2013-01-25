//
//  WellBoreAnalysis.cpp
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 __MyCompanyName__. All rights reserved.
//
#include "WellBoreAnalysis.h"
#include "GeoMeshClass.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic2D.h"
#include "BrazilianTestGeoMesh.h"
#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include <iostream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif
TPZWellBoreAnalysis::TPZWellBoreAnalysis() : fCurrentConfig(), fPostProcessNumber(0)
{
    
}

TPZWellBoreAnalysis::~TPZWellBoreAnalysis()
{
    
}

void TPZWellBoreAnalysis::StandardConfiguration(TPZWellBoreAnalysis &obj)
{
    TPZGeoMesh *gmesh = &obj.fCurrentConfig.fGMesh;
    GeoMeshClass::WellBore2d(&obj.fCurrentConfig.fGMesh);
    obj.fCurrentConfig.fInnerRadius = 0.1;
    obj.fCurrentConfig.fOuterRadius = 1.;
    ofstream arg("wellgeomeshlog.txt");
    obj.fCurrentConfig.fGMesh.Print(arg);
    
    obj.fCurrentConfig.fCMesh.SetReference(gmesh);
    obj.fCurrentConfig.fCMesh.SetDefaultOrder(2);
    TPZCompMesh *compmesh1 = &obj.fCurrentConfig.fCMesh;
    
    TPZSandlerDimaggio::UncDeepSandTest(obj.fCurrentConfig.fSD);
    TPZSandlerDimaggio &SD = obj.fCurrentConfig.fSD;
    SD.fResTol = 1.e-11;
    SD.fIntegrTol = 1.;
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
	TPZMatElastoPlastic2D<TPZSandlerDimaggio> *PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio>(1,1);
    
    obj.fCurrentConfig.fConfinement.XX() = -44.3;
    obj.fCurrentConfig.fConfinement.YY() = -58.2;
    obj.fCurrentConfig.fConfinement.ZZ() = -53.8;
    
    
    TPZTensor<REAL> initstress,finalstress;
    REAL hydro = obj.fCurrentConfig.fConfinement.I1();
    hydro -= SD.fYC.fA*SD.fYC.fR;
    hydro /= 3.;
    finalstress.XX() = hydro;
    finalstress.YY() = hydro;
    finalstress.ZZ() = hydro;
    
    PrepareInitialMat(SD, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = obj.fCurrentConfig.fConfinement;
    PrepareInitialMat(SD, initstress, finalstress, 10);
    
    TPZMaterial *plastic(PlasticSD);
    
    PlasticSD->SetPlasticity(SD);
    compmesh1->InsertMaterialObject(plastic);
    
    CmeshWell(compmesh1,plastic,obj.fCurrentConfig.fConfinement);

}

TPZWellBoreAnalysis::TConfig::TConfig() : fInnerRadius(0.), fOuterRadius(0.), fConfinement(), fSD(), fFluidPressure(0.),
                    fGMesh(), fCMesh(), fAllSol(), fPlasticDeformSqJ2()
{
}

TPZWellBoreAnalysis::TConfig::TConfig(const TConfig &conf) : fInnerRadius(conf.fInnerRadius), fOuterRadius(conf.fOuterRadius), 
        fConfinement(conf.fConfinement), fSD(conf.fSD), fFluidPressure(conf.fFluidPressure),
        fGMesh(conf.fGMesh), fCMesh(conf.fCMesh), fAllSol(conf.fAllSol), fPlasticDeformSqJ2(conf.fPlasticDeformSqJ2)
{
    fGMesh.ResetReference();
    fCMesh.SetReference(&fGMesh);
    fCMesh.LoadReferences();
}

TPZWellBoreAnalysis::TConfig::~TConfig()
{
    
}

TPZWellBoreAnalysis::TConfig &TPZWellBoreAnalysis::TConfig::operator=(const TPZWellBoreAnalysis::TConfig &copy)
{
    fInnerRadius = copy.fInnerRadius;
    fOuterRadius = copy.fOuterRadius;
    fConfinement = copy.fConfinement;
    fSD = copy.fSD;
    fFluidPressure = copy.fFluidPressure;
    fGMesh = copy.fGMesh;
    fCMesh = copy.fCMesh;
    fCMesh.SetReference(&fGMesh);
    fAllSol = copy.fAllSol;
    fPlasticDeformSqJ2 = copy.fPlasticDeformSqJ2;
    return *this;
}

void TPZWellBoreAnalysis::ExecuteInitialSimulation()
{
    TPZElastoPlasticAnalysis analysis(&fCurrentConfig.fCMesh,std::cout);

	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(4);
	analysis.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
	analysis.SetSolver(step);

    int BCId=-2;
    TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    int nsteps = 5;
    
#ifdef LOG4CXX
    {
        TPZMaterial *matmain = analysis.Mesh()->FindMaterial(1);
        std::stringstream sout;
        matmain->Print(sout);
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
        if (!pMatWithMem) {
            DebugStop();
        }

        pMatWithMem->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    TPZFNMatrix<9,STATE> mattemp(3,3,0.), matinit(3,3,0.),matfinal(3,3,0.), matincrement(3,3,0.);
    matinit = pBC->Val1();
    matfinal(0,0) = -29.3;
    matfinal(1,1) = -29.3;
    
    matincrement = matfinal-matinit;
    matincrement *= (1./nsteps);
    mattemp = matinit;
    
    
        
    int neq = analysis.Mesh()->Solution().Rows();
    fCurrentConfig.fAllSol.Redim(neq, nsteps+1);
    
    
    for(int i=0;i<=nsteps;i++)
    {
        pBC->Val1()=mattemp;
        analysis.IterativeProcess(cout, 1.e-8, 30);
        
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
        
        
        
        TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            fCurrentConfig.fAllSol(ieq,i) = sol(ieq,0);
        }
        
        analysis.AcceptSolution();
        
        fCurrentConfig.ComputeElementDeformation();
        
        if (i==0) {
            PostProcess(0);// pOrder
        }
        else {
            PostProcess(1);// pOrder
        }
        
        fPostProcessNumber++;            
        
        mattemp += matincrement;
        
    }
    
    fSequence.push_back(fCurrentConfig);
#ifdef LOG4CXX
    std::list<TConfig>::reverse_iterator it = fSequence.rbegin();
    std::stringstream sout;
    it->fCMesh.Print(sout);
    LOGPZ_DEBUG(logger, sout.str())
#endif
}

void TPZWellBoreAnalysis::ExecuteSimulation()
{
    TPZElastoPlasticAnalysis analysis(&fCurrentConfig.fCMesh,std::cout);
    
	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(4);
	analysis.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
    step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    
    
#ifdef LOG4CXX
    {
        TPZMaterial *matmain = analysis.Mesh()->FindMaterial(1);
        std::stringstream sout;
        matmain->Print(sout);
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
        if (!pMatWithMem) {
            DebugStop();
        }
        
        pMatWithMem->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    
    int neq = analysis.Mesh()->Solution().Rows();
    fCurrentConfig.fAllSol.Redim(neq, 1);
    
    
    analysis.IterativeProcess(cout, 1.e-8, 30);
        
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
        
        
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        fCurrentConfig.fAllSol(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();
    
    fCurrentConfig.ComputeElementDeformation();
    
    PostProcess(1);// pOrder
    fPostProcessNumber++;

    fSequence.push_back(fCurrentConfig);
}


/// Apply the deformation of the configuration to the element
void TPZWellBoreAnalysis::TConfig::ApplyDeformation(TPZCompEl *cel)
{
    TPZCompEl *cel2 = cel;
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    if (!intel2) {
        DebugStop();
    }
    TPZCompMesh *cmesh2 = cel2->Mesh();
    TPZGeoMesh *gmesh1 = &fGMesh;
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh2->MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }
    
    if (intel2->Material() != pMatWithMem2) {
        DebugStop();
    }
    const TPZIntPoints &intpoints = intel2->GetIntegrationRule();
    int nint = intpoints.NPoints();
    TPZManVector<REAL,3> point(3,0.);
    TPZMaterialData data1,data2;
    intel2->InitMaterialData(data2);
    data2.fNeedsSol = false;
    data1.fNeedsSol = true;
    int data2phir = data2.phi.Rows();
    int data2phic = data2.phi.Cols();
    int data2dphir = data2.dphix.Rows();
    int data2dphic = data2.dphix.Cols();
    int elementid = 0;
    TPZManVector<REAL,3> qsi(2,0.);

    for (int ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data2.intLocPtIndex = ip;
        intel2->ComputeRequiredData(data2, point);
        if (data2.intGlobPtIndex == 8334) {
            std::cout << "I should stop\n";
        }
#ifdef LOG4CXX
        {
            int memoryindex = data2.intGlobPtIndex;
            std::stringstream sout;
            sout << "Local integration point index " << data2.intLocPtIndex << std::endl;
            pMatWithMem2->PrintMem(sout,memoryindex);
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif
        TPZGeoEl *gel1 = gmesh1->FindElement(data2.x, qsi, elementid,2);
        if (!gel1) {
            DebugStop();
        }
        TPZCompEl *cel1 = gel1->Reference();
        if (!cel1) {
            DebugStop();
        }
        TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
        if (!intel1) {
            DebugStop();
        }
        intel1->InitMaterialData(data1);
        data1.fNeedsSol = true;
        intel1->ComputeRequiredData(data1, qsi);
#ifdef DEBUG
        {
            REAL diff = dist(data1.x,data2.x);
            if(diff > 1.e-6)
            {
                DebugStop();
            }
        }
#endif
        data2.sol = data1.sol;
        data2.dsol = data1.dsol;
        data2.phi.Redim(0, 1);
        data2.dphix.Redim(2, 0);
        TPZFMatrix<STATE> ek,ef;
        pMatWithMem2->Contribute(data2, weight, ek, ef);
        data2.phi.Resize(data2phir, data2phic);
        data2.dphix.Resize(data2dphir, data2dphic);
    }
    pMatWithMem2->SetUpdateMem(false);

}

void TPZWellBoreAnalysis::TransferSolutionTo(TConfig &config)
{
    TPZCompMesh &cmesh2 = config.fCMesh;
    TPZCompMesh &cmesh1 = fCurrentConfig.fCMesh;
    
    cmesh1.Solution() = fCurrentConfig.fAllSol;
    
    int nel = cmesh2.NElements();
    TPZGeoMesh *gmesh1 = cmesh1.Reference();
    TPZMaterial *mat1 = cmesh1.FindMaterial(1);
    if (!mat1) {
        DebugStop();
    }
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh2.MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }

    
    int varindex1 = mat1->VariableIndex("Stress");
    if (varindex1 == -1) {
        DebugStop();
    }
    int elementid = 0;
    TPZManVector<REAL,3> qsi(3,0.);
    
    for (int el =0; el<nel; el++) {
        TPZCompEl *cel2 = cmesh2.ElementVec()[el];
        if (!cel2) {
            continue;
        }
        TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
        if (!intel2) {
            continue;
        }
        if (intel2->Material() != pMatWithMem2) {
            continue;
        }
        fCurrentConfig.ApplyDeformation(intel2);
    }
    pMatWithMem2->SetUpdateMem(false);
}

// Get the vector of element plastic deformations
void TPZWellBoreAnalysis::TConfig::ComputeElementDeformation()
{
    int nelem = fCMesh.NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fCMesh.ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
        DebugStop();
    }
    for (int el = 0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        fPlasticDeformSqJ2[el] = 0.;
        if (!cel) {
            continue;
        }
        TPZManVector<int> memindices;
        cel->GetMemoryIndices(memindices);
        int numind = memindices.size();
        REAL sqj2el = 0.;
        for (int ind=0; ind<numind; ind++) 
        {
            int memoryindex = memindices[ind];
            if (memoryindex < 0) {
                continue;
            }
            TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
            TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
            REAL J2 = plastic.J2();
            REAL sqj2 = sqrt(J2);
            sqj2el = max(sqj2,sqj2el);
        }
        fPlasticDeformSqJ2[el] = sqj2el;
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Element deformation " << fPlasticDeformSqJ2;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    fCMesh.SetElementSolution(0, fPlasticDeformSqJ2);
}

void TPZWellBoreAnalysis::TConfig::DeleteElementsAbove(REAL sqj2)
{
    int nelem = fCMesh.NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        if (fPlasticDeformSqJ2[el] > sqj2) {
            delete cel;
        }
    }
    // put boundary conditions on the sides which have no neighbours
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel) {
            DebugStop();
        }
        int ns = gel->NSides();
        for (int is=0; is<ns; is++) {
            if (gel->SideDimension(is) != 1) {
                continue;
            }
            TPZStack<TPZCompElSide> cneighbours;
            TPZGeoElSide gelside(gel,is);
            gelside.ConnectedCompElementList(cneighbours, 1, 0);
            if (cneighbours.size() == 0) {
                if (gel->Dimension() == 1) {
                    delete cel;
                }
                else if (gel->Dimension() == 2) {
                    // create a boundary condition
                    TPZGeoElBC gbc(gelside, -3);
                    TPZGeoEl *bc = gbc.CreatedElement();
                    // create the corresponding computational element
                    int index;
                    fCMesh.CreateCompEl(bc, index);
                }
            }
        }
    }
}

/// Change the polynomial order of element using the plastic deformation as threshold
void TPZWellBoreAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<int> &elindices)
{
    DebugStop();
}

/// Divide the element using the plastic deformation as threshold
void TPZWellBoreAnalysis::TConfig::DivideElementsAbove(REAL sqj2, std::set<int> &elindices)
{
    fGMesh.ResetReference();
    fCMesh.LoadReferences();
    int nelem = fCMesh.NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        if (!intel) {
            DebugStop();
        }
        if (fCMesh.ElementSolution()(el,0) < sqj2) {
            continue;
        }
        int porder = intel->GetPreferredOrder();
        TPZStack<int> subels;
        int index = cel->Index();
        intel->Divide(index, subels,0);
        for (int is=0; is<subels.size(); is++) {
            elindices.insert(subels[is]);
            TPZCompEl *subcel = fCMesh.ElementVec()[subels[is]];
            TPZInterpolationSpace *subintel = dynamic_cast<TPZInterpolationSpace *>(subcel);
            if (!subintel) {
                DebugStop();
            }
            subintel->SetPreferredOrder(porder);
        }
    }
    fCMesh.AdjustBoundaryElements();
    fCMesh.InitializeBlock();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// Reset the plastic memory of the integration points of these elements
void TPZWellBoreAnalysis::ApplyHistory(std::set<int> &elindices)
{
    std::list<TConfig>::iterator listit;
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        listit->fCMesh.LoadSolution(listit->fAllSol);
    }

    std::set<int>::iterator it;
    for (it=elindices.begin(); it != elindices.end(); it++) {
        int elindex = *it;
        TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[elindex];
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            DebugStop();
        }
        // Reset the memory of the integration points of the element
        TPZManVector<int> pointindices;
        cel->GetMemoryIndices(pointindices);
        int npoints = pointindices.size();
        for (int ip = 0; ip<npoints; ip++) {
            int ind = pointindices[ip];
            pMatWithMem2->ResetMemItem(ind);
        }
        std::list<TConfig>::iterator listit;
        for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
            listit->ApplyDeformation(cel);
        }
    }
}


void TPZWellBoreAnalysis::PostProcess(int resolution)
{
    std::string vtkFile = "pocoplastico.vtk";
    TPZPostProcAnalysis ppanalysis(&fCurrentConfig.fCMesh);
    ppanalysis.SetStep(fPostProcessNumber);
    TPZFStructMatrix structmatrix(ppanalysis.Mesh());
    structmatrix.SetNumThreads(8);
    ppanalysis.SetStructuralMatrix(structmatrix);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZStack<std::string> PostProcVars, scalNames, vecNames;
    scalNames.Push("Alpha");
    scalNames.Push("PlasticSqJ2");
    scalNames.Push("PlasticSqJ2El");
    //    scalnames[1] = "PlasticSteps";
    //    scalnames[2] = "VolElasticStrain";
    //    scalnames[3] = "VolPlasticStrain";
    //    scalnames[4] = "VolTotalStrain";
    scalNames.Push("I1Stress");
    scalNames.Push("J2Stress");
    
	//vecNames.Push("Displacement");
    vecNames.Push("YieldSurface");
    vecNames.Push("NormalStress");
    vecNames.Push("ShearStress");
    vecNames.Push("NormalStrain");
    vecNames.Push("ShearStrain");
    for (int i=0; i<scalNames.size(); i++) {
        PostProcVars.Push(scalNames[i]);
    }
    for (int i=0; i<vecNames.size(); i++) {
        PostProcVars.Push(vecNames[i]);
    }
    //    
    ppanalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
    //
    ppanalysis.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
    //	
    ppanalysis.TransferSolution();
    ppanalysis.PostProcess(resolution);// pOrder
    fPostProcessNumber++;

}

/// Divide the element using the plastic deformation as threshold
void TPZWellBoreAnalysis::DivideElementsAbove(REAL sqj2)
{
    std::set<int> elindices;
    fCurrentConfig.DivideElementsAbove(sqj2,elindices);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Element indices that will be refined ";
        for (std::set<int>::iterator it = elindices.begin(); it!= elindices.end(); it++) {
            sout << *it << " ";
        }
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fCurrentConfig.fCMesh.Print(sout);    
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // subject the integration points with the deformation history
    ApplyHistory(elindices);
    fCurrentConfig.ComputeElementDeformation();
}
