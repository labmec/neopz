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
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>

#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif


void CmeshWell(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure)
{
/*    TPZFMatrix<REAL> BeginStress(3,3,0.), EndStress(2,2,0.), EndStress2(3,3,0.);
    TPZFMatrix<REAL> val1(3,1,0.);TPZFMatrix<REAL> val2(3,1,0.);TPZFMatrix<REAL> BeginForce(3,1,0.);TPZFMatrix<REAL> EndForce(3,1,0.);
    TPZTensor<REAL> OCStress, beginOCStress, loadStress, loadStress2, initialStrain, FarFieldStress, TestStress;
	
	REAL alpha=0.8;
	REAL PorePressure = 4397.;//PSI
    REAL L=3000.,LDA=1200.,terra=0.9,agua=0.447;//PSI/ft
    REAL SigmaV  = agua*3.28*LDA+terra*3.28*(L-LDA);
	REAL Sigmah = 0.8 * (SigmaV - PorePressure * alpha) + PorePressure * alpha;
    
	loadStress.XX()=1.;
    loadStress.YY()=1.;
	
	loadStress *= 1.*(Sigmah - PorePressure * alpha);
    EndStress(0,0)=1.*(Sigmah - PorePressure * alpha);
    EndStress(1,1)=1.*(Sigmah - PorePressure * alpha);
    loadStress.CopyToTensor(EndStress);

 */
    
    
    TPZFMatrix<REAL> f1(3,1,0.);
    TPZFMatrix<REAL> k1(3,3,0.);
    k1(0,0)=Confinement.XX();
    k1(1,1)=Confinement.YY();
    k1(2,2)=Confinement.ZZ();
    TPZMaterial *bc1 = mat->CreateBC(mat,-2,4,k1,f1);
    CMesh->InsertMaterialObject(bc1);

    
    // type 6 constraints in x and y
    // type 5 only normal constraint
    TPZFNMatrix<9> k5(3,3,0.),f5(3,1,pressure);
    for (int i=0; i<3; i++) {
        k5(i,i) = 1.e5;
    }
    TPZMaterial *bc5 = mat->CreateBC(mat, -6, 6, k5, f5);
    CMesh->InsertMaterialObject(bc5);

   
    TPZFMatrix<REAL> k2(3,3,0.);
    TPZFMatrix<REAL> f2(3,1,0.);
    k2(0,0)=Confinement.XX();
    k2(1,1)=Confinement.YY();
    k2(2,2)=Confinement.ZZ();
    TPZMaterial * bc2 = mat->CreateBC(mat,-3,4,k2,f2);
    CMesh->InsertMaterialObject(bc2);
    
    
    TPZFMatrix<REAL> k3(2,2,0.);
    TPZFMatrix<REAL> f3(2,1,0.);
    f3(1,0)=1.;
   // k3(1,1)=1.e12;
    TPZMaterial * bc3 = mat->CreateBC(mat,-4,3,k3,f3);
    CMesh->InsertMaterialObject(bc3);
    
    TPZFMatrix<REAL> k4(2,2,0.);
    TPZFMatrix<REAL> f4(2,1,0.);
     f4(0,0)=1.;
    //k4(0,0)=1.e12;
    
    TPZMaterial * bc4 = mat->CreateBC(mat,-5,3,k4,f4);
    CMesh->InsertMaterialObject(bc4);
//    CMesh->SetDefaultOrder(2);
    CMesh->AutoBuild();

#ifdef LOG4CXX
    {
        std::stringstream sout;
        CMesh->Print(sout);
        LOGPZ_DEBUG(logger , sout.str())
    }
#endif
}


TPZWellBoreAnalysis::TPZWellBoreAnalysis() : fCurrentConfig(), fSequence(), fPostProcessNumber(0)
{
    
}

TPZWellBoreAnalysis::~TPZWellBoreAnalysis()
{
    
}

/// write the object on the stream
void TPZWellBoreAnalysis::Write(TPZStream &output)
{
    fCurrentConfig.Write(output);
    int seqsize = fSequence.size();
    output.Write(&seqsize);
    std::list<TPZWellBoreAnalysis::TConfig>::iterator it;
    for (it = fSequence.begin(); it != fSequence.end(); it++) {
        it->Write(output);
    }
    output.Write(&fPostProcessNumber);
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fCurrentConfig.fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
}

/// read the object from the stream
void TPZWellBoreAnalysis::Read(TPZStream &input)
{
    fCurrentConfig.Read(input);
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fCurrentConfig.fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    int seqsize;
    input.Read(&seqsize);
    for (int i=0; i<seqsize; i++) {
        TPZWellBoreAnalysis::TConfig config;
        config.Read(input);
        fSequence.push_back(config);
    }
    input.Read(&fPostProcessNumber);
    
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
    obj.fCurrentConfig.fCMesh.SetDefaultOrder(1);
    TPZCompMesh *compmesh1 = &obj.fCurrentConfig.fCMesh;
    
//    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::UncDeepSandTest(obj.fCurrentConfig.fSD);
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::PRSMatMPa(obj.fCurrentConfig.fSD);
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &SD = obj.fCurrentConfig.fSD;
    SD.SetResidualTolerance(1.e-10);
    SD.fIntegrTol = 10.;
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
	TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > *PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> >(1,1);
    
    obj.fCurrentConfig.fConfinement.XX() = -44.3;// MPa
    obj.fCurrentConfig.fConfinement.YY() = -58.2;
    obj.fCurrentConfig.fConfinement.ZZ() = -53.8;
    obj.fCurrentConfig.fFluidPressure = -29.3;
    
    
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
    SD.ResetPlasticMem();
    
    TPZMaterial *plastic(PlasticSD);
    
    PlasticSD->SetPlasticity(SD);
    compmesh1->InsertMaterialObject(plastic);
    
    CmeshWell(compmesh1,plastic,obj.fCurrentConfig.fConfinement,obj.fCurrentConfig.fFluidPressure);

}

void TPZWellBoreAnalysis::CheckDeformation(std::string filename)
{
    
    std::ofstream out(filename.c_str());
    TPZTensor<STATE> epstotal,sigma, S, pressure ;
    sigma.XX() = 1.;
    sigma.XY() = 1.;
    STATE i1 = sigma.I1();
    pressure.Identity();
    pressure *= i1/3.;
    sigma.S(S);
    STATE sqj2 = sqrt(sigma.J2());
    S *= 1./sqj2;
    sqj2 = S.J2();
    const int nangles = 100;
    out << "deform = Table[0,{" << nangles << "}];" << std::endl;
    out << "tension = Table[0,{" << nangles << "}];" << std::endl;
    
    for (int i=0; i<nangles; i++) {
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> SD;
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::PRSMatMPa(SD);
        SD.SetResidualTolerance(1.e-10);
        SD.fIntegrTol = 1.;

        const int nincrements = 100;
        TPZManVector<std::pair<STATE, STATE>,100 > epst(nincrements),sigvec(nincrements);
        STATE angle = i*M_PI/nangles;
        STATE ca = cos(angle);
        STATE sa = sin(angle);
        TPZTensor<STATE> deformbase,p1(pressure), s1(S);
        p1 *= ca;
        s1 *= sa;
        deformbase = p1;
        deformbase += s1;
        STATE catest = deformbase.I1();
        STATE satest = sqrt(deformbase.J2());
        catest -= ca;
        satest -= fabs(sa);
        for (int j=0; j< nincrements; j++) {
            STATE scale = 0.1*j/nincrements;
            TPZTensor<STATE> epstotal(deformbase);
            epstotal *= scale;
            TPZTensor<STATE> sigma;
            STATE i1 = epstotal.I1();
            STATE sqj2deform = sqrt(epstotal.J2());
            STATE diff[2] = {i1-ca*scale,sqj2deform-fabs(sa)*scale};
            if (fabs(diff[0]) > 1.e-6 || fabs(diff[1]) > 1.e-6) {
                DebugStop();
            }
            epst[j].first = ca*scale;
            epst[j].second = sa*scale;
            SD.fYC.fIsonCap = false;
            SD.ApplyStrainComputeSigma(epstotal,sigma);
            // print I1 and sqrt(J2)
            REAL i1sigma = sigma.I1();
            REAL sqj2 = sqrt(sigma.J2());
            sigvec[j].first = i1sigma;
            sigvec[j].second = sqj2;
        }
        out << "deform[[" << i+1 << "]] = " <<  epst << ";" << std::endl;
        out << "tension[[" << i+1 << "]] = " << sigvec << ";" << std::endl;
    }

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
    if (fGMesh.NElements() > 87) 
    {
        std::cout << __PRETTY_FUNCTION__ << std::endl;
        fGMesh.ElementVec()[87]->Print();
    }
    fCMesh.CleanUp();
    fCMesh.SetReference(0);
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

/// Write the data to the output stream
void TPZWellBoreAnalysis::TConfig::Write(TPZStream &out)
{
    out.Write(&fInnerRadius);
    out.Write(&fOuterRadius);
    fConfinement.Write(out);
    fSD.Write(out);
    out.Write(&fFluidPressure);
    fGMesh.Write(out, 0);
    fCMesh.Write(out, 0);
    fAllSol.Write(out, 0);
    TPZSaveable::WriteObjects(out,fPlasticDeformSqJ2);
    int verify = 83562;
    out.Write(&verify);
}

/// Read the data from the input stream
void TPZWellBoreAnalysis::TConfig::Read(TPZStream &input)
{
    input.Read(&fInnerRadius);
    input.Read(&fOuterRadius);
    fConfinement.Read(input);
    fSD.Read(input);
    input.Read(&fFluidPressure);
    fGMesh.Read(input, 0);
    fCMesh.Read(input, &fGMesh);
    fAllSol.Read(input, 0);
    TPZSaveable::ReadObjects(input,fPlasticDeformSqJ2);
    int verify = 0;
    input.Read(&verify);
    if (verify != 83562)
    {
        DebugStop();
    }
}



void TPZWellBoreAnalysis::ExecuteInitialSimulation()
{
    TPZElastoPlasticAnalysis analysis(&fCurrentConfig.fCMesh,std::cout);

	TPZSkylineStructMatrix full(&fCurrentConfig.fCMesh);
    full.SetNumThreads(4);
	analysis.SetStructuralMatrix(full);
    
    analysis.AddNoPenetration(-5, 0);
    analysis.AddNoPenetration(-4, 1);
    analysis.IdentifyEquationsToZero();
    
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
    matfinal(0,0) = fCurrentConfig.fFluidPressure;//-29.3;
    matfinal(1,1) = fCurrentConfig.fFluidPressure;//-29.3;
    
    matincrement = matfinal-matinit;
    matincrement *= (1./nsteps);
    mattemp = matinit;
    
    
        
    int neq = analysis.Mesh()->Solution().Rows();
    fCurrentConfig.fAllSol.Redim(neq, nsteps+1);
    
    
    for(int i=0;i<=nsteps;i++)
    {
        std::cout << "Initial Simulation Step " << i << " out of " << nsteps << std::endl;
        pBC->Val1()=mattemp;
        int numNewton = 7;
//        bool linesearch = true;
//        bool checkconv = true;
//        analysis.IterativeProcess(cout, 1.e-8, numNewton,linesearch,checkconv);
        analysis.IterativeProcess(cout, 1.e-8, numNewton);
        
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
        
        
        
        TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            fCurrentConfig.fAllSol(ieq,i) = sol(ieq,0);
        }
        
        TPZPlasticDiagnostic diag(analysis.Mesh());
        diag.CheckGlobal();
        
        analysis.AcceptSolution();
        
        fCurrentConfig.ComputeElementDeformation();
        
        if (i==0) {
            PostProcess(0);// pOrder
        }
        else {
            PostProcess(1);// pOrder
        }
        
        fCurrentConfig.VerifyGlobalEquilibrium();
        //fPostProcessNumber++;            
        
        mattemp += matincrement;
        
    }
    
    fCurrentConfig.VerifyGlobalEquilibrium();
    
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
    
    analysis.AddNoPenetration(-5, 0);
    analysis.AddNoPenetration(-4, 1);
    analysis.IdentifyEquationsToZero();

	TPZStepSolver<REAL> step;
    step.SetDirect(ELDLt);
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
    
    int NumIter = 7;
    analysis.IterativeProcess(cout, 1.e-8, NumIter);
        
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
    
    fCurrentConfig.VerifyGlobalEquilibrium();

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

/// Verify the global equilibrium of the forces by boundary condition
void TPZWellBoreAnalysis::TConfig::VerifyGlobalEquilibrium(std::ostream &out)
{
    int neq = fCMesh.NEquations();
    
    TPZFMatrix<STATE> rhsDirichlet(neq,1,0.);
    this->ComputeRhsExceptMatid(-4, rhsDirichlet);
//    this->ComputeRhsExceptMatid(-5, rhsDirichlet);
    
    TPZStack<int,10> boundaryconditionlist;
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        TPZMaterial *mat = itmat->second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        int matid = itmat->first;
        allmaterials.insert(matid);
        if (matid == -4 || matid == -5) {
            continue;
        }
        if (bnd) {
            boundaryconditionlist.Push(itmat->first);
        }
    }
    std::map<int, TPZManVector<REAL,2> > ForceResultantsExcept, ForceResultantsOnly;
    for(int ib=0; ib<boundaryconditionlist.size(); ib++)
    {
        int bc = boundaryconditionlist[ib];
        {
            TPZFMatrix<STATE> rhs(neq,1,0.);
            this->ComputeRhsForMatid(bc, rhs);
            this->ComputeXYForce(rhs , ForceResultantsOnly[ib]);
        }
        {
            TPZFMatrix<STATE> rhs(neq,1,0.);
            this->ComputeRhsExceptMatid(bc, rhs);
            rhs -= rhsDirichlet;
            this->ComputeXYForce(rhs , ForceResultantsExcept[ib]);
        }
        
    }
    
    TPZManVector<REAL,2> totalforce(2,0.), dirichletforce(2,0.);
    ComputeXYForce(rhsDirichlet, dirichletforce);
    totalforce[0] = -dirichletforce[0];
    totalforce[1] = -dirichletforce[1];
    out <<  "FORCE RESULTANTS\n";
    out << "Dirichlet force " << dirichletforce << std::endl;
    for (int ib=0; ib<boundaryconditionlist.size(); ib++) {
        int bc = boundaryconditionlist[ib];
        if (bc<0) {
            totalforce[0] += ForceResultantsOnly[ib][0];
            totalforce[1] += ForceResultantsOnly[ib][1];
        }
        out << "Force computing only BC = " << bc << " x-force = " << ForceResultantsOnly[ib][0] << " y-force = " << ForceResultantsOnly[ib][1] << std::endl;
        out << "Force computing except BC = " << bc << " x-force = " << ForceResultantsExcept[ib][0] << " y-force = " << ForceResultantsExcept[ib][1] << std::endl;
    }
    out << "Force resultants " << totalforce << std::endl;
}

/// Verify the global equilibrium of the forces by boundary condition
void TPZWellBoreAnalysis::TConfig::VerifyGlobalEquilibrium2(std::ostream &out)
{
    TPZStack<int,10> boundaryconditionlist;
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        TPZMaterial *mat = itmat->second;
        TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
        allmaterials.insert(itmat->first);
        if (bnd) {
            boundaryconditionlist.Push(itmat->first);
        }
    }
    if (boundaryconditionlist[0] != -6 || boundaryconditionlist[2] != -4) {
        DebugStop();
    }
    boundaryconditionlist[0] = -4;
    boundaryconditionlist[2] = -6;
    TPZManVector<TPZManVector<REAL,2>, 10> ForceResultants(boundaryconditionlist.size()+1);
    int neq = fCMesh.NEquations();
    TPZFMatrix<STATE> rhs4(neq,1,0.),rhs5(neq,1,0.);
    for(int ib=0; ib<boundaryconditionlist.size(); ib++)
    {
        int bc = boundaryconditionlist[ib];
        TPZFStructMatrix str(&fCMesh);
        std::set<int> matset(allmaterials);
        matset.erase(bc);
        str.SetMaterialIds(matset);
        TPZFMatrix<STATE> rhs(neq,1,0.);
        //        if (bc > -4) {
        //            matset.erase(-5);
        //            matset.erase(-4);
        //        }
        str.Assemble(rhs,0);
        // correct for the Dirichlet boundary conditions
        rhs -= rhs4;
        rhs -= rhs5;
        //        if (bc == -5) {
        //            rhs5 = rhs;
        //        }
        //        if (bc == -4) {
        //            rhs4 = rhs;
        //        }
        std::set<int> xeqs,yeqs;
        for (int el=0; el<fCMesh.NElements(); el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel || !cel->Reference()) {
                continue;
            }
            int matid = cel->Reference()->MaterialId();
            if (matid != bc) {
                continue;
            }
            int nc = cel->Reference()->NCornerNodes();
            for (int ic=0; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                int seq = c.SequenceNumber();
                int pos = fCMesh.Block().Position(seq);
                xeqs.insert(pos);
                yeqs.insert(pos+1);
            }
        }
        ForceResultants[ib].Resize(2,0.);
        std::set<int>::iterator it;
        for(it = xeqs.begin(); it != xeqs.end(); it++)
        {
            ForceResultants[ib][0] += rhs(*it,0);
            if (bc == -5) {
                rhs5(*it,0) = rhs(*it,0);
            }
        }
        for(it = yeqs.begin(); it != yeqs.end(); it++)
        {
            ForceResultants[ib][1] += rhs(*it,0);
            if (bc == -4) {
                rhs4(*it,0) = rhs(*it,0);
            }
        }
    }
    REAL normrhs = 0.;
    {
        boundaryconditionlist.Push(1);
        int ib = boundaryconditionlist.size()-1;
        ForceResultants[ib].Resize(2,0.);
        TPZFStructMatrix str(&fCMesh);
        TPZFMatrix<STATE> rhs(neq,1,0.);
        str.Assemble(rhs, 0);
        rhs -= rhs4;
        rhs -= rhs5;
        normrhs = Norm(rhs);
        ForceResultants[ib].Resize(2, 0.);
        for (int i=0; i<neq/2; i++) {
            ForceResultants[ib][0] += rhs(2*i,0);
            ForceResultants[ib][1] += rhs(2*i+1,0);
        }
    }
    
    TPZManVector<REAL,2> totalforce(2,0.);
    out <<  "FORCE RESULTANTS\n";
    for (int ib=0; ib<boundaryconditionlist.size(); ib++) {
        int bc = boundaryconditionlist[ib];
        if (bc<0) {
            totalforce[0] += ForceResultants[ib][0];
            totalforce[1] += ForceResultants[ib][1];
        }
        out << "Boundary condition " << bc << " x-force = " << ForceResultants[ib][0] << " y-force = " << ForceResultants[ib][1] << std::endl;
    }
    out << "Force resultants " << totalforce << std::endl;
    out << "Norm rhs " << normrhs << std::endl;
}

/// Compute the Rhs for the mesh minus the elements with matid
// this method is cumulative (sums to the rhs)
void TPZWellBoreAnalysis::TConfig::ComputeRhsExceptMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    std::map<int,TPZMaterial *>::iterator itmat;
    std::map<int,TPZMaterial *> &matmap = fCMesh.MaterialVec();
    for (itmat = matmap.begin(); itmat != matmap.end(); itmat++) {
        allmaterials.insert(itmat->first);
    }
    allmaterials.erase(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetMaterialIds(allmaterials);
    str.Assemble(rhs, 0);

    
}

/// Compute the contribution of only matid
// this method is cumulative (sums to the rhs)
void TPZWellBoreAnalysis::TConfig::ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    allmaterials.insert(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetMaterialIds(allmaterials);
    str.Assemble(rhs, 0);
}

/// Compute the resultant x and y force
void TPZWellBoreAnalysis::TConfig::ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force)
{
    force.Resize(2, 0.);
    force.Fill(0.);
    int nel = rhs.Rows();
    for (int i=0; i<nel; i++) {
        force[i%2] += rhs(i,0);
    }
}


void TPZWellBoreAnalysis::TConfig::DeleteElementsAbove(REAL sqj2)
{
    fCMesh.Reference()->ResetReference();
    fCMesh.LoadReferences();
    int ndel = 0;
    int nelem = fCMesh.NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        if (fPlasticDeformSqJ2[el] > sqj2) {
            TPZGeoEl *gel = cel->Reference();
            delete cel;
            gel->SetMaterialId(50);
            ndel++;
        }
        //gel->ResetReference();
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
                    TPZGeoElBC gbc(gelside, -6);
                    TPZGeoEl *bc = gbc.CreatedElement();
                    // create the corresponding computational element
                    int index;
                    fCMesh.CreateCompEl(bc, index);
                }
            }
        }
    }
    TPZMaterial *mat = fCMesh.FindMaterial(-6);
    if (!mat) {
        DebugStop();
    }
    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
    if (!bnd) {
        DebugStop();
    }
    bnd->Val1()(0,0) = 1.e8;
    bnd->Val1()(1,1) = 1.e8;
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    {
        std::ofstream file("AdjustedMesh.vtk");
        TPZVTKGeoMesh::PrintGMeshVTK(fCMesh.Reference(), file, true);
    }
    
    std::cout << "Number of elements deleted " << ndel << std::endl;
}

// Diminish the spring which stabilizes the well by the given factor
void TPZWellBoreAnalysis::TConfig::RelaxWellSpring(REAL factor)
{
    TPZMaterial *mat = fCMesh.FindMaterial(-6);
    if (!mat) {
        DebugStop();
    }
    TPZBndCond *bnd = dynamic_cast<TPZBndCond *>(mat);
    if (!bnd) {
        DebugStop();
    }
    bnd->Val1() *= factor;
    bnd->Val1().Print("Value of the spring",std::cout);
    bnd->Val2().Print("Fluid pressure ",std::cout);
}


/// Change the polynomial order of element using the plastic deformation as threshold
void TPZWellBoreAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<int> &elindices)
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
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            continue;
        }

        if (fCMesh.ElementSolution()(el,0) < sqj2) {
            continue;
        }
        TPZStack<int> subels;
        int index = cel->Index();
        elindices.insert(index);
        intel->SetPreferredOrder(porder);
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

TPZGeoMesh * TPZWellBoreAnalysis::TConfig::GetGeoMesh() {
  return &this->fGMesh;
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
    scalNames.Push("POrder");
    //    scalnames[1] = "PlasticSteps";
    //    scalnames[2] = "VolElasticStrain";
    //    scalnames[3] = "VolPlasticStrain";
    //    scalnames[4] = "VolTotalStrain";
    scalNames.Push("I1Stress");
    scalNames.Push("J2Stress");
    
	vecNames.Push("Displacement");
    vecNames.Push("YieldSurface");
    vecNames.Push("NormalStress");
    vecNames.Push("ShearStress");
    vecNames.Push("NormalStrain");
    vecNames.Push("ShearStrain");
    vecNames.Push("DisplacementMem");
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

/// GetPostProcessedValues
void TPZWellBoreAnalysis::PostProcessedValues(TPZVec<REAL> &x, TPZVec<std::string> &variables, TPZFMatrix<STATE> &values)
{
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
    pMatWithMem2->SetUpdateMem(true);
    
    // identify the post process indices and number of results for each index
    int varsize = variables.size();
    TPZManVector<int,10> varindices(varsize,-1),firstindex(varsize+1,-1),numvar(varsize,-1);
    int numvalid = 0;
    firstindex[0] = 0;
    for (int iv=0; iv<varsize; iv++) {
        varindices[numvalid] = pMatWithMem2->VariableIndex(variables[iv]);
        if (varindices[numvalid] != -1) {
            numvar[numvalid] = pMatWithMem2->NSolutionVariables(varindices[numvalid]);
            firstindex[numvalid+1] = firstindex[numvalid]+numvar[numvalid];
            numvalid++;
        }
    }
    varindices.Resize(numvalid);
    firstindex.Resize(numvalid+1);
    numvar.Resize(numvalid);
    values.Resize(0, firstindex[numvalid]);
    
    // compute the number of post processing variables and their first index
    std::list<TConfig>::iterator listit;
    // create a new integration point to train the deformation history
    int pointid = pMatWithMem2->PushMemItem();
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        TPZGeoMesh *gmesh1 = &(listit->fGMesh);
        TPZManVector<REAL,3> qsi(2,0.);
        int elementid;
        TPZMaterialData data1;
        data1.x = x;
        data1.intLocPtIndex = 0;
        // find the element which contains the coordinate
        
        TPZGeoEl *gel1 = gmesh1->FindElement(data1.x, qsi, elementid,2);
        if (!gel1) {
            DebugStop();
        }
        TPZCompEl *cel1 = gel1->Reference();
        if (!cel1) {
            DebugStop();
        }
        // compute the solutions at the point
        TPZInterpolationSpace *intel1 = dynamic_cast<TPZInterpolationSpace *>(cel1);
        if (!intel1) {
            DebugStop();
        }

        // expand the result matrix
        int valuesfirstindex = values.Rows();
        int nsol = listit->fAllSol.Cols();
        // we don't initialize the values because they will be set in the procedure
        values.Resize(valuesfirstindex+nsol, values.Cols());
        // copy the solution to a local variable
        int solsize = listit->fAllSol.Rows();
        TPZFMatrix<STATE> locsol(solsize,1);
        for (int isol = 0; isol < nsol; isol++) 
        {
            // load one solution at a time
            for (int ic=0; ic<solsize; ic++) {
                locsol(ic,0) = listit->fAllSol(ic,isol);
            }
            listit->fCMesh.LoadSolution(locsol);
            intel1->InitMaterialData(data1);
            data1.fNeedsSol = true;
            intel1->ComputeRequiredData(data1, qsi);
            data1.intGlobPtIndex = pointid;
            // initialize the material data object
            TPZFMatrix<STATE> EF(0,0);
            data1.phi.Resize(0, 1);
            data1.dphix.Resize(data1.dphix.Rows(), 0);
            pMatWithMem2->Contribute(data1, 1., EF);
            // call the post processing method
            for (int iv=0; iv<numvalid; iv++) 
            {
                TPZManVector<STATE,10> solution(numvar[iv],0.);
                pMatWithMem2->Solution(data1, varindices[iv], solution);
                // put the data in the output matrix
                for (int is=0; is<numvar[iv]; is++) {
                    values(valuesfirstindex+isol,firstindex[iv]+is) = solution[is];
                }
            }
        }
    }
    pMatWithMem2->FreeMemItem(pointid);
    pMatWithMem2->SetUpdateMem(false);
    
}

