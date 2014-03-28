//
//  WellBoreAnalysis.cpp
//  PZ
//
//  Created by phil on 1/18/13.
//  Copyright (c) 2013 LabMeC. All rights reserved.
//
#include "WellBoreAnalysis.h"
#include "GeoMeshClass.h"
#include "pzelastoplasticanalysis.h"
#include "pzelastoplastic2D.h"

#include "BrazilianTestGeoMesh.h"
#include "tpzchangeel.h"
//#include "poroelastoplastic.h"
#include "pzgeoelbc.h"
#include "TPZVTKGeoMesh.h"
#include "pzplasticdiagnostic.h"
#include <iostream>
#include "pzbfilestream.h"

#include "pzlog.h"


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.plasticity.wellboreanalysis"));
#endif


void CmeshWell(TPZCompMesh *CMesh, TPZMaterial * mat, TPZTensor<STATE> &Confinement, STATE pressure)
{
    
    TPZMaterial *prev = 0;
    
    TPZFMatrix<REAL> f1(3,1,0.);
    TPZFMatrix<REAL> k1(3,3,0.);
    k1(0,0)=Confinement.XX();
    k1(1,1)=Confinement.YY();
    k1(2,2)=Confinement.ZZ();
    TPZMaterial *bc1 = mat->CreateBC(mat,-2,4,k1,f1);
    prev = CMesh->FindMaterial(-2);
    if(prev)
    {
        delete prev;
    }
    CMesh->InsertMaterialObject(bc1);

    
    // type 6 constraints in x and y
    // type 5 only normal constraint
    TPZFNMatrix<9> k5(3,3,0.),f5(3,1,pressure);
    for (int i=0; i<3; i++) {
        k5(i,i) = 1.e5;
    }
    TPZMaterial *bc5 = mat->CreateBC(mat, -6, 6, k5, f5);
    prev = CMesh->FindMaterial(-6);
    if(prev)
    {
        delete prev;
    }
    CMesh->InsertMaterialObject(bc5);

   
    TPZFMatrix<REAL> k2(3,3,0.);
    TPZFMatrix<REAL> f2(3,1,0.);
    k2(0,0)=Confinement.XX();
    k2(1,1)=Confinement.YY();
    k2(2,2)=Confinement.ZZ();
    TPZMaterial * bc2 = mat->CreateBC(mat,-3,4,k2,f2);
    
    prev = CMesh->FindMaterial(-3);
    if(prev)
    {
        delete prev;
    }
 
    CMesh->InsertMaterialObject(bc2);
    
    
    TPZFMatrix<REAL> k3(2,2,0.);
    TPZFMatrix<REAL> f3(2,1,0.);
    f3(1,0)=1.;
   // k3(1,1)=1.e12;
    TPZMaterial * bc3 = mat->CreateBC(mat,-4,3,k3,f3);
    prev = CMesh->FindMaterial(-4);
    if(prev)
    {
        delete prev;
    }

    CMesh->InsertMaterialObject(bc3);
    
    TPZFMatrix<REAL> k4(2,2,0.);
    TPZFMatrix<REAL> f4(2,1,0.);
     f4(0,0)=1.;
    //k4(0,0)=1.e12;
    
    TPZMaterial * bc4 = mat->CreateBC(mat,-5,3,k4,f4);
    prev = CMesh->FindMaterial(-5);
    if(prev)
    {
        delete prev;
    }

    CMesh->InsertMaterialObject(bc4);
//    CMesh->SetDefaultOrder(2);
//    CMesh->AutoBuild();

#ifdef LOG4CXX
    if(logger->isDebugEnabled())
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
    if(logger->isDebugEnabled())
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
    if(logger->isDebugEnabled())
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
    //GeoMeshClass::WellBore2d(&obj.fCurrentConfig.fGMesh);
    obj.fCurrentConfig.fInnerRadius = 4.25*0.0254;//0.1;
    obj.fCurrentConfig.fOuterRadius = 3.;//1.;
    ofstream arg("wellgeomeshlog.txt");
    obj.fCurrentConfig.fGMesh.Print(arg);
    
    obj.fCurrentConfig.fCMesh.SetReference(gmesh);
    int defaultporder = 2;
    obj.fCurrentConfig.fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &obj.fCurrentConfig.fCMesh;
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
#ifdef PV

    int materialid=1;
    bool planestrain=true;
    //obj.fCurrentConfig.fSDPV.fYC.PreSMat(obj.fCurrentConfig.fSDPV.fYC);
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> &SD =obj.fCurrentConfig.fSDPV;
    REAL poisson = 0.203;
    REAL elast = 29269.;
    REAL A = 152.54;
    REAL B = 0.0015489;
    REAL C = 146.29;
    REAL R = 0.91969;
    REAL D = 0.018768;
    REAL W = 0.006605;
    REAL G=elast/(2.*(1.+poisson));
    REAL K=elast/(3.*(1.-2*poisson));
    REAL phi=0,psi=1.,N=0.;
    SD.fYC.SetUp( A,  B, C,  D, K, G, W, R, phi, N, psi);
    SD.fER.SetUp(elast,poisson);
    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> > *PlasticSD = new TPZMatElastoPlastic2D< TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> >(materialid,planestrain);

#else
    
    //TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::PRSMatMPa(obj.fCurrentConfig.fSD);
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &SD = obj.fCurrentConfig.fSD;
    SD.SetResidualTolerance(1.e-10);
    SD.fIntegrTol = 10.;

   REAL poisson = 0.203;
   REAL elast = 29269.;
   REAL A = 152.54;
   REAL B = 0.0015489;
   REAL C = 146.29;
   REAL R = 0.91969;
   REAL D = 0.018768;
   REAL W = 0.006605;
   SD.SetUp(poisson, elast, A, B, C, R, D, W);

	TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > *PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> >(1,1);
    
#endif

    
    obj.fCurrentConfig.fConfinement.XX() = -45.9;//-44.3;// MPa
    obj.fCurrentConfig.fConfinement.YY() = -62.1;//-58.2;
    obj.fCurrentConfig.fConfinement.ZZ() = -48.2;//-53.8;
    obj.fCurrentConfig.fFluidPressure = 19.5;//29.3;
    

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
    compmesh1->AutoBuild();
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
        
#ifdef PV
        TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > SD;
        SD.fYC.PreSMat(SD.fYC);
#else
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> SD;
        TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::PRSMatMPa(SD);
        SD.SetResidualTolerance(1.e-10);
        SD.fIntegrTol = 1.;
#endif


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



TPZWellBoreAnalysis::TConfig::TConfig() : fInnerRadius(0.), fOuterRadius(0.), fNx(2,0),fDelx(0.),fGreater(),fSmaller(),fConfinement(), fSD(), fFluidPressure(0.),
    fGMesh(), fCMesh(), fAllSol(), fPlasticDeformSqJ2(), fHistoryLog()
#ifdef PV
  , fSDPV()
#endif
{
}

TPZWellBoreAnalysis::TConfig::TConfig(const TConfig &conf) : fInnerRadius(conf.fInnerRadius), fOuterRadius(conf.fOuterRadius),fNx(conf.fNx),fDelx(conf.fDelx),
    fGreater(conf.fGreater),fSmaller(conf.fSmaller),
        fConfinement(conf.fConfinement), fSD(conf.fSD), fFluidPressure(conf.fFluidPressure),
        fGMesh(conf.fGMesh), fCMesh(conf.fCMesh), fAllSol(conf.fAllSol), fPlasticDeformSqJ2(conf.fPlasticDeformSqJ2), fHistoryLog(conf.fHistoryLog)
{
    fGMesh.ResetReference();
    fCMesh.SetReference(&fGMesh);
    fCMesh.LoadReferences();
}

TPZWellBoreAnalysis::TConfig::~TConfig()
{
    fPostprocess.SetCompMesh(0);
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
    fPostprocess.SetCompMesh(0);
    fGMesh = copy.fGMesh;
    fCMesh = copy.fCMesh;
    fCMesh.SetReference(&fGMesh);
    fAllSol = copy.fAllSol;
    fPlasticDeformSqJ2 = copy.fPlasticDeformSqJ2;
    fHistoryLog = copy.fHistoryLog;
    return *this;
}

/// Write the data to the output stream
void TPZWellBoreAnalysis::TConfig::Write(TPZStream &out)
{
    out.Write(&fInnerRadius);
    out.Write(&fOuterRadius);
    out.Write(&fDelx);
    TPZSaveable::WriteObjects(out, fNx);
    TPZSaveable::WriteObjects(out, fGreater);
    TPZSaveable::WriteObjects(out, fSmaller);
    fConfinement.Write(out);
    fSD.Write(out);
    out.Write(&fFluidPressure);
    fGMesh.Write(out, 0);
    fCMesh.Write(out, 0);
    fAllSol.Write(out, 0);
    TPZSaveable::WriteObjects(out,fPlasticDeformSqJ2);
    out.Write(&fHistoryLog);

    int verify = 83562;
    out.Write(&verify);
}

/// Read the data from the input stream
void TPZWellBoreAnalysis::TConfig::Read(TPZStream &input)
{
    input.Read(&fInnerRadius);
    input.Read(&fOuterRadius);
    input.Read(&fDelx);
    TPZSaveable::ReadObjects(input, fNx);
    TPZSaveable::ReadObjects(input, fGreater);
    TPZSaveable::ReadObjects(input, fSmaller);
    
    fConfinement.Read(input);
    fSD.Read(input);
    input.Read(&fFluidPressure);
    fGMesh.Read(input, 0);
    fCMesh.Read(input, &fGMesh);
    fAllSol.Read(input, 0);
    TPZSaveable::ReadObjects(input,fPlasticDeformSqJ2);
    input.Read(&fHistoryLog);
    
    int verify = 0;
    input.Read(&verify);
    if (verify != 83562)
    {
        DebugStop();
    }
}



void TPZWellBoreAnalysis::ExecuteInitialSimulation(int nsteps, int numnewton)
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
    //step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    
//    analysis.Assemble();

    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
    if (!pMatWithMem) {
        DebugStop();
    }
    pMatWithMem->ResetMemory();
    
    CmeshWell(&fCurrentConfig.fCMesh, pMatWithMem, fCurrentConfig.fConfinement, fCurrentConfig.fFluidPressure);
    
    int BCId=-2;
    TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);

    TPZFNMatrix<9,STATE> mattemp(3,3,0.), matinit(3,3,0.),matfinal(3,3,0.), matincrement(3,3,0.);
    matinit = pBC->Val1();
    matfinal(0,0) = -fCurrentConfig.fFluidPressure;//-29.3;
    matfinal(1,1) = -fCurrentConfig.fFluidPressure;//-29.3;
    
    matincrement = matfinal-matinit;
    matincrement *= (1./nsteps);
    mattemp = matinit;
    
    
        
    int neq = analysis.Mesh()->Solution().Rows();
    fCurrentConfig.fAllSol.Redim(neq, nsteps+1);
    
    
    for(int i=0;i<=nsteps;i++)
    {
        std::cout << "Initial Simulation Step " << i << " out of " << nsteps << std::endl;
        pBC->Val1()=mattemp;
        
        if (fLinearMatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << fLinearMatrix->IsDecomposed() << std::endl;
        }
        
        
        ComputeLinearMatrix();        
        analysis.AdjustTangentMatrix(fLinearMatrix);
        analysis.Solver().SetMatrix(fLinearMatrix);
        if (fLinearMatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << fLinearMatrix->IsDecomposed() << std::endl;
        }
        bool linesearch = true,checkconv=false;
		REAL tol = 1.e-6;
#ifdef PV
        //analysis.IterativeProcess(cout, 1.e-8, numNewton,linesearch,checkconv);
        analysis.IterativeProcess(cout, fLinearMatrix,tol, numnewton, linesearch);
#else
        //analysis.IterativeProcess(cout, 1.e-8, numNewton,linesearch,checkconv);
        analysis.IterativeProcess(cout, fLinearMatrix, tol, numnewton, linesearch);
#endif
        
        

        
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
        if (fLinearMatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << fLinearMatrix->IsDecomposed() << std::endl;
        }

        
        
        TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            fCurrentConfig.fAllSol(ieq,i) = sol(ieq,0);
        }
        
//        TPZPlasticDiagnostic diag(analysis.Mesh());
//        diag.CheckGlobal();
//        
//        if (i == nsteps) {
//            VerifyTangentValidity();
//        }
        
        analysis.AcceptSolution();
        
        if (fLinearMatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << fLinearMatrix->IsDecomposed() << std::endl;
        }

        
        fCurrentConfig.ComputeElementDeformation();
        

        
        if (i==0) {
            fCurrentConfig.CreatePostProcessingMesh(fPostProcessNumber, 0);
            PostProcess();
        }
        else {
            fCurrentConfig.CreatePostProcessingMesh(fPostProcessNumber,1);
            PostProcess();
        }
        
        if (fLinearMatrix) {
            std::cout << __FILE__ << ":" << __LINE__ << "Decomposed " << fLinearMatrix->IsDecomposed() << std::endl;
        }

        //fCurrentConfig.VerifyGlobalEquilibrium();
        //fPostProcessNumber++;            

        std::stringstream strout;
        strout << "Step " << fPostProcessNumber-1;
        fCurrentConfig.fHistoryLog = strout.str();

        fSequence.push_back(fCurrentConfig);

        mattemp += matincrement;
        
    }
    
//    fCurrentConfig.VerifyGlobalEquilibrium();
    
//    fSequence.push_back(fCurrentConfig);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) 
    {
        std::list<TConfig>::reverse_iterator it = fSequence.rbegin();
        std::stringstream sout;
        it->fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
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
    //step.SetDirect(ECholesky);
	analysis.SetSolver(step);
    
    ComputeLinearMatrix();
    analysis.AdjustTangentMatrix(fLinearMatrix);
    analysis.Solver().SetMatrix(fLinearMatrix);
    
    
    
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
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
    
    int NumIter = 60;
    bool linesearch = true;
    bool checkconv = false;
//    analysis.IterativeProcess(cout, 1.e-6, NumIter, linesearch, checkconv);
    analysis.IterativeProcess(cout, fLinearMatrix, 1.e-8, NumIter, linesearch);
    
    
        //analysis.Solution().Print();
        //analysis.AcceptSolution();
        //analysis.TransferSolution(ppanalysis);
        
//    VerifyTangentValidity();    
    
        
    TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
    for (int ieq=0; ieq<neq; ieq++) {
        fCurrentConfig.fAllSol(ieq,0) = sol(ieq,0);
    }
    
    analysis.AcceptSolution();
    
    fCurrentConfig.ComputeElementDeformation();
    
    fCurrentConfig.CreatePostProcessingMesh(fPostProcessNumber, 1);
    
    PostProcess();
    
//    fCurrentConfig.VerifyGlobalEquilibrium();

    std::stringstream strout;
    strout << "Step " << fPostProcessNumber-1 << " pwb=" << fCurrentConfig.fFluidPressure;
    fCurrentConfig.fHistoryLog = strout.str();

    fSequence.push_back(fCurrentConfig);
}



void TPZWellBoreAnalysis::ExecuteSimulation(int nsteps,REAL pwb)
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
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCurrentConfig.fCMesh.MaterialVec()[1]);
    if (!pMatWithMem) {
        DebugStop();
    }
//    pMatWithMem->ResetMemory();
    
//    CmeshWell(&fCurrentConfig.fCMesh, pMatWithMem, fCurrentConfig.fConfinement, fCurrentConfig.fFluidPressure);
    
    int BCId=-2;
    TPZMaterial * mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    
    TPZFNMatrix<9,STATE> mattemp(3,3,0.), matinit(3,3,0.),matfinal(3,3,0.), matincrement(3,3,0.);
    matfinal(0,0) = -pwb;
    matfinal(1,1) = -pwb;
    matinit(0,0) = -fCurrentConfig.fFluidPressure;//-29.3;
    matinit(1,1) = -fCurrentConfig.fFluidPressure;//-29.3;
    matinit(2,2)=0.;
    
    matincrement = matfinal-matinit;
    matincrement *= (1./nsteps);
    mattemp = matinit + matincrement;
    
    int neq = analysis.Mesh()->Solution().Rows();
    fCurrentConfig.fAllSol.Redim(neq, nsteps+1);
    
    
    for(int i=1;i<=nsteps;i++)
    {
        std::cout << "Simulation Step " << i << " out of " << nsteps << std::endl;
        pBC->Val1()=mattemp;
        
        ComputeLinearMatrix();
        analysis.AdjustTangentMatrix(fLinearMatrix);
        analysis.Solver().SetMatrix(fLinearMatrix);
        
        bool linesearch = true;
        bool checkconv=false;
        int numNewton =60;
#ifdef PV
        //analysis.IterativeProcess(cout, 1.e-8, numNewton,linesearch,checkconv);
        analysis.IterativeProcess(cout, fLinearMatrix, 1.e-8, numNewton, linesearch);
#else

        analysis.IterativeProcess(cout, fLinearMatrix, 1.e-8, numNewton, linesearch);
#endif
        
        TPZFMatrix<STATE> &sol = analysis.Mesh()->Solution();
        for (int ieq=0; ieq<neq; ieq++) {
            fCurrentConfig.fAllSol(ieq,i) = sol(ieq,0);
        }
        
//        TPZPlasticDiagnostic diag(analysis.Mesh());
//        diag.CheckGlobal();
//        
//        if (i == nsteps) {
//            VerifyTangentValidity();
//        }
        analysis.AcceptSolution();
        fCurrentConfig.ComputeElementDeformation();

        
        if (i==0) {
            fCurrentConfig.CreatePostProcessingMesh(fPostProcessNumber, 0);
            PostProcess();
        }
        else {
            fCurrentConfig.CreatePostProcessingMesh(fPostProcessNumber, 1);
            PostProcess();
        }
        
        cout << "-------------------> i: "<< i << " Pressao atual: " << mattemp;

        //fCurrentConfig.VerifyGlobalEquilibrium();

        fCurrentConfig.fFluidPressure=-mattemp(0,0);
        //fCurrentConfig.VerifyGlobalEquilibrium();

        std::stringstream strout;
        strout << "Step " << fPostProcessNumber-1 << " pwb=" << fCurrentConfig.fFluidPressure;
        fCurrentConfig.fHistoryLog = strout.str();

        fSequence.push_back(fCurrentConfig);

        mattemp += matincrement;
    }

//    fCurrentConfig.fFluidPressure=pwb;
//    fCurrentConfig.VerifyGlobalEquilibrium();
//    fSequence.push_back(fCurrentConfig);

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
    TPZCompMesh *cmesh1 = gmesh1->Reference();
    
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh2->MaterialVec()[1]);
    if (pMatWithMem2) {
        pMatWithMem2->SetUpdateMem(true);
    }
    else {
        DebugStop();
    }
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem1 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cmesh1->MaterialVec()[1]);
    
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
    long data2phir = data2.phi.Rows();
    long data2phic = data2.phi.Cols();
    long data2dphir = data2.dphix.Rows();
    long data2dphic = data2.dphix.Cols();
    long elementid = 0;
    TPZManVector<REAL,3> qsi(2,0.);

    for (long ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data2.intLocPtIndex = ip;
        intel2->ComputeRequiredData(data2, point);
#ifdef LOG4CXX2
        if(logger->isDebugEnabled())
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
#ifdef LOG4CXX
        if(logger->isDebugEnabled())
        {
            int memoryindex = data2.intGlobPtIndex;
            std::stringstream sout;
            sout << "Local integration point index " << data2.intLocPtIndex << std::endl;
            sout << "qsi coordinate " << qsi << std::endl;
            pMatWithMem2->PrintMem(sout,memoryindex);
            sout << "\noriginal element index " << cel2->Index() << std::endl;
            sout << "projected element index " << cel1->Index() << std::endl;
            if(cel1->Index() == cel2->Index())
            {
                sout << "Same point index in the original mesh " << std::endl;
                sout << "qsi coordinate " << qsi << std::endl;
                pMatWithMem1->PrintMem(sout,memoryindex);
            }
            LOGPZ_DEBUG(logger, sout.str())
        }
#endif

    }
    pMatWithMem2->SetUpdateMem(false);

}

/// verify if the stress tangent is computed correctly
void TPZWellBoreAnalysis::VerifyTangentValidity()
{
    TPZCompMesh &cmesh = fCurrentConfig.fCMesh;
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        cmesh.Solution().Print("Mesh solution ",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    int nelem = cmesh.NElements();
    for (int el=0; el<nelem; el++) {
        TPZCompEl *cel = cmesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if (!gel || gel->MaterialId() != 1) {
            continue;
        }
        fCurrentConfig.VerifyPlasticTangent(cel);
    }
}


/// Verify tangent elasto plastic relation
void TPZWellBoreAnalysis::TConfig::VerifyPlasticTangent(TPZCompEl *cel)
{
    TPZCompEl *cel2 = cel;
    TPZInterpolationSpace *intel2 = dynamic_cast<TPZInterpolationSpace *>(cel2);
    if (!intel2) {
        DebugStop();
    }
    TPZCompMesh *cmesh2 = cel2->Mesh();
    
    
#ifdef PV
    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *pMatWithMem2 =
    dynamic_cast<TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse > > *>(cmesh2->MaterialVec()[1]);
#else
    TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > *pMatWithMem2 =
    dynamic_cast<TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > *>(cmesh2->MaterialVec()[1]);
#endif

    
    if (intel2->Material() != pMatWithMem2) {
        DebugStop();
    }
    const TPZIntPoints &intpoints = intel2->GetIntegrationRule();
    int nint = intpoints.NPoints();
    TPZManVector<REAL,3> point(3,0.);
    TPZMaterialData data2;
    intel2->InitMaterialData(data2);
    data2.fNeedsSol = true;
    TPZManVector<REAL,3> qsi(2,0.);
    std::stringstream sout;
    sout << "Diagnostic for element " << cel->Index() << std::endl;
    
    
    for (int ip =0; ip<nint; ip++) {
        REAL weight;
        intpoints.Point(ip, point, weight);
        data2.intLocPtIndex = ip;
        intel2->ComputeRequiredData(data2, point);
        TPZFNMatrix<36,REAL> deltastrain(3,1),stress(3,1),dstressdstrain(3,3);
        pMatWithMem2->ComputeDeltaStrainVector(data2, deltastrain);
        pMatWithMem2->ApplyDeltaStrainComputeDep(data2, deltastrain, stress, dstressdstrain);
        TPZFMatrix<REAL> variation(deltastrain);
        variation *= -1/10.;
        TPZManVector<STATE,20> errors(9),convrate(8,0.);
        REAL minerror = 0.;
        for (int ist=1; ist<10; ist++) {
            TPZFNMatrix<6,REAL> deltastrainnext(3,1),stressnext(3,1),stressestimate(3,1),stresserror(3,1);
            for (int i=0; i<3; i++) deltastrainnext(i,0) = deltastrain(i,0)+ist*variation(i,0);
            pMatWithMem2->ApplyDeltaStrain(data2, deltastrainnext, stressnext);
            for (int i=0; i<3; i++) {
                stressestimate(i,0) = stress(i,0);
                for (int j=0; j<3; j++) {
                    stressestimate(i,0) += dstressdstrain(i,j)*ist*variation(j,0);
                }
            }
            for (int i=0; i<3; i++) stresserror(i,0) = stressnext(i,0)-stressestimate(i,0);
            errors[ist-1] = Norm(stresserror);
            if(ist == 1) minerror = errors[ist-1];
            if (minerror > errors[ist-1]) {
                minerror = errors[ist-1];
            }
        }
        if(minerror > 1.e-12)
        {
            for (int ist=2; ist<10; ist++) {
                convrate[ist-2] = log(errors[ist-1]/errors[0])/log(ist*1.);
            }
        }
        sout << "Integration point " << ip << " coordinate " << data2.x << std::endl;
        sout << "Errors " << errors << std::endl;
        sout << "Convergence rates " << convrate << std::endl;
    }
#ifdef LOG4CXX
    LOGPZ_DEBUG(logger, sout.str())
#else
    std::cout << sout.str();
#endif
}



void TPZWellBoreAnalysis::TransferSolutionTo(TConfig &config)
{
    TPZCompMesh &cmesh2 = config.fCMesh;
    TPZCompMesh &cmesh1 = fCurrentConfig.fCMesh;
    
    cmesh1.Solution() = fCurrentConfig.fAllSol;
    
    long nel = cmesh2.NElements();
//    TPZGeoMesh *gmesh1 = cmesh1.Reference();
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
//    int elementid = 0;
    TPZManVector<REAL,3> qsi(3,0.);
    
    for (long el =0; el<nel; el++) {
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

/// Compute the area of the domain at which sqJ2 is above a given value
REAL TPZWellBoreAnalysis::TConfig::ComputeAreaAboveSqJ2(REAL sqj2)
{
    REAL area = 0.;
    long nelem = fCMesh.NElements();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel) {
                continue;
            }
            bool shouldanalyse = false;
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            long numind = memindices.size();
  //          REAL sqj2el = 0.;
            for (long ind=0; ind<numind; ind++) 
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2pt = sqrt(J2);
                if (sqj2pt > sqj2) {
                    shouldanalyse = true;
                }
            }
            if (shouldanalyse) {
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                TPZGeoEl *gel = cel->Reference();
                if (! intel || !gel) {
                    DebugStop();
                }
                TPZIntPoints &rule = intel->GetIntegrationRule();
                int np = rule.NPoints();
                for (int ip = 0; ip<np; ip++) {
                    TPZManVector<REAL,3> point(2,0.);
                    REAL weight;
                    rule.Point(ip, point, weight);
                    TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
                    TPZFNMatrix<9,REAL> axes(2,3);
                    REAL detjac;
                    gel->Jacobian(point, jac, axes, detjac, jacinv);
                    TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ip]);
                    TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                    REAL J2 = plastic.J2();
                    REAL sqj2pt = sqrt(J2);
                    if (sqj2pt > sqj2) {
                        area += weight*fabs(detjac);
                    }
                }
            }
        }
    }
    return area;
}

/// Compute the area of the domain at which sqJ2 is above a given value
REAL TPZWellBoreAnalysis::TConfig::OpeningAngle(REAL sqj2)
{
    REAL angle = 0.;
    long nelem = fCMesh.NElements();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel) {
                continue;
            }
            bool shouldanalyse = false;
            TPZManVector<long> memindices;
            cel->GetMemoryIndices(memindices);
            int numind = memindices.size();
            REAL sqj2el;
            sqj2el = 0.;
            for (int ind=0; ind<numind; ind++) 
            {
                int memoryindex = memindices[ind];
                if (memoryindex < 0) {
                    continue;
                }
                TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ind]);
                TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                REAL J2 = plastic.J2();
                REAL sqj2pt = sqrt(J2);
                if (sqj2pt > sqj2) {
                    shouldanalyse = true;
                }
            }
            if (shouldanalyse) {
                TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
                TPZGeoEl *gel = cel->Reference();
                if (! intel || !gel) {
                    DebugStop();
                }
                TPZIntPoints &rule = intel->GetIntegrationRule();
                int np = rule.NPoints();
                for (int ip = 0; ip<np; ip++) {
                    TPZManVector<REAL,3> point(2,0.);
                    REAL weight;
                    rule.Point(ip, point, weight);
                    TPZManVector<REAL,3> x(3,0.);
                    gel->X(point, x);
                    TPZElastoPlasticMem &mem = pMatWithMem2->MemItem(memindices[ip]);
                    TPZTensor<REAL> &plastic = mem.fPlasticState.fEpsP;
                    REAL J2 = plastic.J2();
                    REAL sqj2pt = sqrt(J2);
                    if (sqj2pt > sqj2) {
                        REAL locangle = atan2(x[1],x[0]);
                        if (locangle > angle) {
                            angle = locangle;
                        }
                    }
                }
            }
        }
    }
    return angle;
}

/// Compute the removed area of the domain
REAL TPZWellBoreAnalysis::TConfig::RemovedArea()
{
    return -ComputeTotalArea()+M_PI*(fOuterRadius*fOuterRadius-fInnerRadius*fInnerRadius)/4.;
}

/// Compute the area of the domain
REAL TPZWellBoreAnalysis::TConfig::ComputeTotalArea()
{
    REAL area = 0.;
    long nelem = fCMesh.NElements();
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel) {
                continue;
            }
            TPZGeoEl *gel = cel->Reference();
            if (gel->MaterialId() != 1) {
                continue;
            }
            TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
            if (!intel) {
                DebugStop();
            }
            TPZIntPoints &rule = intel->GetIntegrationRule();
            int np = rule.NPoints();
            for (int ip = 0; ip<np; ip++) {
                TPZManVector<REAL,3> point(2,0.);
                REAL weight;
                rule.Point(ip, point, weight);
                TPZFNMatrix<4,REAL> jac(2,2),jacinv(2,2);
                TPZFNMatrix<9,REAL> axes(2,3);
                REAL detjac;
                gel->Jacobian(point, jac, axes, detjac, jacinv);
                area += weight*fabs(detjac);
            }
        }
    }
    return area;    
}


// Get the vector of element plastic deformations
void TPZWellBoreAnalysis::TConfig::ComputeElementDeformation()
{
    long nelem = fCMesh.NElements();
    fPlasticDeformSqJ2.resize(nelem);
    fPlasticDeformSqJ2.Fill(0.);
    fCMesh.ElementSolution().Redim(nelem, 1);
    TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (fCMesh.MaterialVec()[1]);
    if (!pMatWithMem2) {
        fPlasticDeformSqJ2.Fill(0.);
    }
    else 
    {
        for (long el = 0; el<nelem; el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            fPlasticDeformSqJ2[el] = 0.;
            if (!cel) {
                continue;
            }
            TPZManVector<long> memindices;
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
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
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
    long neq = fCMesh.NEquations();
    
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
        std::set<long> xeqs,yeqs;
        for (long el=0; el<fCMesh.NElements(); el++) {
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
                long seq = c.SequenceNumber();
                long pos = fCMesh.Block().Position(seq);
                xeqs.insert(pos);
                yeqs.insert(pos+1);
            }
        }
        ForceResultants[ib].Resize(2,0.);
        std::set<long>::iterator it;
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
        // zerar os residuos referente a connects que nao sao de vertice
        for (long el=0; el<fCMesh.NElements(); el++) {
            TPZCompEl *cel = fCMesh.ElementVec()[el];
            if (!cel || !cel->Reference()) {
                continue;
            }
            int ncorner = cel->Reference()->NCornerNodes();
            long nc = cel->NConnects();
            for (long ic=ncorner; ic<nc; ic++) {
                TPZConnect &c = cel->Connect(ic);
                long seq = c.SequenceNumber();
                long pos = fCMesh.Block().Position(seq);
                int blsize = fCMesh.Block().Size(seq);
                for (int ibl =0; ibl<blsize; ibl++) {
                    rhs(pos+ibl,0) = 0.;
                }
            }
        }
        
        
        rhs -= rhs4;
        rhs -= rhs5;
        normrhs = Norm(rhs);
        ForceResultants[ib].Resize(2, 0.);
        for (long i=0; i<neq/2; i++) {
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

/// Zera os componentes do rhs para connects diferentes do zero
void TPZWellBoreAnalysis::TConfig::FilterRhs(TPZFMatrix<STATE> &rhs)
{
    // zerar os residuos referente a connects que nao sao de vertice
    for (long el=0; el<fCMesh.NElements(); el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel || !cel->Reference()) {
            continue;
        }
        int ncorner = cel->Reference()->NCornerNodes();
        int nc = cel->NConnects();
        for (int ic=ncorner; ic<nc; ic++) {
            TPZConnect &c = cel->Connect(ic);
            if(c.HasDependency())
            {
                continue;
            }
            long seq = c.SequenceNumber();
            long pos = fCMesh.Block().Position(seq);
            int blsize = fCMesh.Block().Size(seq);
            for (int ibl =0; ibl<blsize; ibl++) {
                rhs(pos+ibl,0) = 0.;
            }
        }
    }
    
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
    str.SetNumThreads(8);
    str.SetMaterialIds(allmaterials);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);

    
}

/// Compute the contribution of only matid
// this method is cumulative (sums to the rhs)
void TPZWellBoreAnalysis::TConfig::ComputeRhsForMatid(int matid, TPZFMatrix<STATE> &rhs)
{
    std::set<int> allmaterials;
    allmaterials.insert(matid);
    TPZFStructMatrix str(&fCMesh);
    str.SetMaterialIds(allmaterials);
    str.SetNumThreads(8);
    str.Assemble(rhs, 0);
    FilterRhs(rhs);
}

/// Compute the resultant x and y force
void TPZWellBoreAnalysis::TConfig::ComputeXYForce(TPZFMatrix<STATE> &rhs, TPZVec<STATE> &force)
{
    force.Resize(2, 0.);
    force.Fill(0.);
    long nel = rhs.Rows();
    for (long i=0; i<nel; i++) {
        force[i%2] += rhs(i,0);
    }
}


void TPZWellBoreAnalysis::TConfig::DeleteElementsAbove(REAL sqj2)
{
    fCMesh.Reference()->ResetReference();
    fCMesh.LoadReferences();
    long ndel = 0;
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
        TPZCompEl *cel = fCMesh.ElementVec()[el];
        if (!cel) {
            continue;
        }
        if (fPlasticDeformSqJ2[el] > sqj2) {
            TPZGeoEl *gel = cel->Reference();
            // delete the neighbouring boundary elements
            int ncorner = gel->NCornerNodes();
            int ns = gel->NSides();
            if(gel->Index() == 1212)
            {
                std::cout << __PRETTY_FUNCTION__ << " I should stop\n";
                gel->Print(std::cout);
            }
            for (int is=ncorner; is<ns; is++) {
                TPZGeoElSide neighbour = gel->Neighbour(is);
                TPZGeoEl *neighgel = neighbour.Element();
                TPZCompEl *neighcel = neighgel->Reference();
                long neighcelindex = 0;
                if(neighcel) neighcelindex = neighcel->Index();
                if (neighgel->Dimension() == 1 && neighcel) {
                    delete neighcel;
                    neighgel->SetMaterialId(50);
                }
                if(gel->Index() == 1212)
                {
                    std::cout << __PRETTY_FUNCTION__ << " I should stop\n";
                    neighgel->Print(std::cout);
                    std::cout << fCMesh.ElementVec()[neighcelindex] << std::endl;
                }
            }
            delete cel;
            gel->SetMaterialId(50);
            ndel++;
        }
        //gel->ResetReference();
    }
    // put boundary conditions on the sides which have no neighbours
    for (long el=0; el<nelem; el++) {
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
                    gel->SetMaterialId(50);
                }
                else if (gel->Dimension() == 2) {
                    // create a boundary condition
                    TPZGeoElBC gbc(gelside, -6);
                    TPZGeoEl *bc = gbc.CreatedElement();
                    // create the corresponding computational element
                    long index;
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
void TPZWellBoreAnalysis::TConfig::PRefineElementsAbove(REAL sqj2, int porder, std::set<long> &elindices)
{
    fGMesh.ResetReference();
    fCMesh.LoadReferences();
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
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
        TPZStack<long> subels;
        long index = cel->Index();
        elindices.insert(index);
        intel->SetPreferredOrder(porder);
    }
    fCMesh.AdjustBoundaryElements();
    fCMesh.InitializeBlock();
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// Divide the element using the plastic deformation as threshold
void TPZWellBoreAnalysis::TConfig::DivideElementsAbove(REAL sqj2, std::set<long> &elindices)
{
    fGMesh.ResetReference();
    fCMesh.LoadReferences();
    long nelem = fCMesh.NElements();
    for (long el=0; el<nelem; el++) {
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
        int porder = intel->GetPreferredOrder();
        TPZStack<long> subels;
        long index = cel->Index();
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
    if(logger->isDebugEnabled())
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
void TPZWellBoreAnalysis::ApplyHistory(std::set<long> &elindices)
{
    std::list<TConfig>::iterator listit;
    for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
        listit->fCMesh.LoadSolution(listit->fAllSol);
    }
    
    std::stringstream filename;
    filename << "applyhistory_" << startfrom << ".txt";
    std::ofstream out(filename.str().c_str());

    std::set<long>::iterator it;
    for (it=elindices.begin(); it != elindices.end(); it++) {
        long elindex = *it;
        TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[elindex];
        TPZMatWithMem<TPZElastoPlasticMem> *pMatWithMem2 = dynamic_cast<TPZMatWithMem<TPZElastoPlasticMem> *> (cel->Material());
        if (!pMatWithMem2) {
            DebugStop();
        }
        // Reset the memory of the integration points of the element
        TPZManVector<long> pointindices;
        cel->GetMemoryIndices(pointindices);
        long npoints = pointindices.size();
        for (long ip = 0; ip<npoints; ip++) {
            long ind = pointindices[ip];
            pMatWithMem2->ResetMemItem(ind);
        }
        std::list<TConfig>::iterator listit;
        for (listit = fSequence.begin(); listit != fSequence.end(); listit++) {
            listit->ApplyDeformation(cel);
        }
        for (long ip = 0; ip<npoints; ip++) {
            long ind = pointindices[ip];
            pMatWithMem2->MemItem(ind).Print(out);
        }
    }
}

void TPZWellBoreAnalysis::TConfig::CreatePostProcessingMesh(int PostProcessNumber, int resolution)
{
#ifdef PV
    std::string vtkFile = "pocoplasticoPV.vtk";
    //std::string vtkFile = "pocoplasticoPV2.vtk";
#else
    std::string vtkFile = "pocoplasticoErick.vtk";
#endif
    if (fPostprocess.ReferenceCompMesh() == &fCMesh) {
        return;
    }
    
    fPostprocess.SetCompMesh(&fCMesh);
    TPZFStructMatrix structmatrix(fPostprocess.Mesh());
    structmatrix.SetNumThreads(8);
    fPostprocess.SetStructuralMatrix(structmatrix);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZStack<std::string> PostProcVars, scalNames, vecNames;
    scalNames.Push("Alpha");//
    scalNames.Push("PlasticSqJ2");//
    scalNames.Push("PlasticSqJ2El");//
    scalNames.Push("POrder");//
    //    scalnames[1] = "PlasticSteps";
    //    scalnames[2] = "VolElasticStrain";
    //    scalnames[3] = "VolPlasticStrain";
    //    scalnames[4] = "VolTotalStrain";
    scalNames.Push("I1Stress");//
    scalNames.Push("J2Stress");//
    scalNames.Push("DisplacementX");//
    scalNames.Push("DisplacementY");//
    scalNames.Push("DisplacementZ");//
    scalNames.Push("PrincipalStressDirection1");//
    scalNames.Push("PrincipalStressDirection2");//
    scalNames.Push("PrincipalStressDirection3");//
    scalNames.Push("StressX");//
    scalNames.Push("StressY");//
    scalNames.Push("StressZ");//
#ifdef MustComeBack
    scalNames.Push("I1HorStress");
    scalNames.Push("J2HorStress");
    
    vecNames.Push("Displacement");
    vecNames.Push("YieldSurface");
    vecNames.Push("NormalStress");//
    vecNames.Push("ShearStress");//
    vecNames.Push("NormalStrain");//
    vecNames.Push("ShearStrain");//
    vecNames.Push("DisplacementMem");
    vecNames.Push("PrincipalStress");
#endif
    vecNames.Push("DisplacementMem");//

    for (int i=0; i<scalNames.size(); i++) {
        PostProcVars.Push(scalNames[i]);
    }
    for (int i=0; i<vecNames.size(); i++) {
        PostProcVars.Push(vecNames[i]);
    }
    //
    fPostprocess.SetPostProcessVariables(PostProcMatIds, PostProcVars);
    //
//    fPostprocess.DefineGraphMesh(2,scalNames,vecNames,vtkFile);
    fPostprocess.TransferSolution();
    
    fPostprocess.SetStep(PostProcessNumber);
    if (fPostprocess.ReferenceCompMesh() != &fCMesh) {
        fPostprocess.SetCompMesh(&fCMesh);
    }
//    fPostprocess.PostProcess(resolution);
}

int passCount = 0;
void TPZWellBoreAnalysis::PostProcess()
{
    fPostProcessNumber++;

    //*** DANGER (definido pelo programador) ***
    REAL J2val = 0.0004;
    //******************************************
    
#ifdef TESTPOLYCHAIN
    std::multimap<REAL,REAL> polygonalChain;

    GetJ2Isoline(J2val, polygonalChain);
    
    {
        std::stringstream nm;
        nm << "EllipDots" << passCount << ".nb";
        std::ofstream outEllips(nm.str().c_str());
        passCount++;
        std::multimap<REAL,REAL>::iterator it;
        outEllips << "ellipDots={";
        for(it = polygonalChain.begin(); it != polygonalChain.end(); it++)
        {
            outEllips << "{" << it->first << " , " << it->second << "}";
            if(it != polygonalChain.end())
            {
                outEllips << ",";
            }
        }
        outEllips << "};\n";
        outEllips << "aa=ListPlot[ellipDots,Joined->False,AspectRatio->Automatic];\n";
        outEllips << "bb=Graphics[Circle[{0,0},4.25*0.0254]];\n";
        outEllips << "Show[aa,bb]\n";
    }
#endif
}

void TPZWellBoreAnalysis::GetJ2Isoline(REAL J2val, std::multimap<REAL,REAL> & polygonalChain)
{
    TPZCompMesh *cmesh = fCurrentConfig.fPostprocess.Mesh();
    REAL tol = 0.01*J2val;
    
    long nels = cmesh->NElements();
    for(long el = 0; el < nels; el++)
    {
        TPZCompEl * cel = cmesh->ElementVec()[el];
        if(!cel)
        {
            continue;
        }
        if(cel->Dimension() == 2)
        {
            int nnodes = cel->Reference()->NCornerNodes();
            TPZVec<REAL> nodeSol(nnodes,0.);
            TPZVec< TPZVec<REAL> > qsiNode(nnodes);
            
            TPZVec<STATE> sol(1);
            int var = cel->Material()->VariableIndex("PlasticSqJ2");
            REAL Tol = 1.e-6;
			//ZeroTolerance(Tol);
            for(int n = 0; n < nnodes; n++)
            {
                qsiNode[n].Resize(2, 0.);
                TPZVec<REAL> xNode(3,0.);
                cel->Reference()->NodePtr(n)->GetCoordinates(xNode);
                cel->Reference()->ComputeXInverse(xNode, qsiNode[n],Tol);
                cel->Solution(qsiNode[n], var, sol);
                nodeSol[n] = sol[0];
            }
            for(int n = 0; n < nnodes; n++)
            {
                int pos0 = (n)%nnodes;
                int pos1 = (n+1)%nnodes;
                REAL minSol = MIN(nodeSol[pos0],nodeSol[pos1]);
                REAL maxSol = MAX(nodeSol[pos0],nodeSol[pos1]);
                
                if(J2val - minSol > tol && maxSol - J2val > tol)
                {//Bisseccao!!!
                    TPZVec<REAL> qsiMin(2,0.), qsiMax(2,0.);
                    if(fabs(minSol - nodeSol[pos0]) < tol)
                    {
                        qsiMin = qsiNode[pos0];
                        qsiMax = qsiNode[pos1];
                    }
                    else
                    {
                        qsiMin = qsiNode[pos1];
                        qsiMax = qsiNode[pos0];
                    }
                    TPZVec<REAL> qsiMiddle(2,0.);
                    qsiMiddle[0] = (qsiMin[0] + qsiMax[0])/2.;
                    qsiMiddle[1] = (qsiMin[1] + qsiMax[1])/2.;
                    
                    cel->Solution(qsiMiddle, var, sol);
                    
                    REAL difSol = fabs(sol[0] - J2val);
                    
                    while(difSol > tol)
                    {
                        if(sol[0] < J2val)
                        {
                            qsiMin = qsiMiddle;
                        }
                        else
                        {
                            qsiMax = qsiMiddle;
                        }
                        
                        qsiMiddle[0] = (qsiMin[0] + qsiMax[0])/2.;
                        qsiMiddle[1] = (qsiMin[1] + qsiMax[1])/2.;
                        
                        cel->Solution(qsiMiddle, var, sol);
                        
                        difSol = fabs(sol[0] - J2val);
                    }
                    TPZManVector<REAL,3> searchedX(3,0.);
                    cel->Reference()->X(qsiMiddle, searchedX);
                    polygonalChain.insert( std::make_pair(searchedX[0],searchedX[1]) );
                }
            }
        }
    }
}

/// Modify the geometry of the domain simulating an elliptic breakout
void TPZWellBoreAnalysis::AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis)
{
    fCurrentConfig.AddEllipticBreakout(MaiorAxis, MinorAxis);
    std::set<long> elindices;
    int nel = fCurrentConfig.fCMesh.NElements();
    for (long el=0; el<nel; el++) {
        TPZCompEl *cel = fCurrentConfig.fCMesh.ElementVec()[el];
        TPZGeoEl *gel = 0;
        if(cel) gel = cel->Reference();
        if(cel && gel && gel->MaterialId() == 1)
        {
            elindices.insert(el);
        }
    }
    ApplyHistory(elindices);
}



/// Divide the element using the plastic deformation as threshold
unsigned int TPZWellBoreAnalysis::DivideElementsAbove(REAL sqj2)
{
    std::set<long> elindices;
    fCurrentConfig.DivideElementsAbove(sqj2,elindices);
#ifdef LOG4CXX
    if (logger->isDebugEnabled()) 
    {
        std::stringstream sout;
        sout << "Element indices that have been created ";
        for (std::set<long>::iterator it = elindices.begin(); it!= elindices.end(); it++) {
            sout << *it << " ";
        }
        std::cout << sout.str() << std::endl;
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fGMesh.Print(sout);    
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fCurrentConfig.fCMesh.Print(sout);    
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif

    // subject the integration points with the deformation history
    ApplyHistory(elindices);
    fCurrentConfig.ComputeElementDeformation();
		return elindices.size();
}

/// GetPostProcessedValues at a given coordinate x
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
        long elementid = 0;
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

/// Add elliptic breakout
void TPZWellBoreAnalysis::TConfig::AddEllipticBreakout(REAL MaiorAxis, REAL MinorAxis)
{
    if (fGreater.size() && fGreater[fGreater.size()-1] > MaiorAxis) {
        DebugStop();
    }
    if (MaiorAxis < fInnerRadius) {
        DebugStop();
    }
    fGreater.Resize(fGreater.size()+1, MaiorAxis);
    if (fSmaller.size() && fSmaller[fSmaller.size()-1] < MinorAxis) {
        DebugStop();
    }
    if (MinorAxis > fInnerRadius) {
        DebugStop();
    }
    fSmaller.Resize(fSmaller.size()+1, MinorAxis);
    
    fPostprocess.SetCompMesh(0);
    fCMesh.CleanUp();
    fGMesh.CleanUp();
    CreateMesh();

    CreateComputationalMesh(2);
    int BCId=-2;
    TPZMaterial * mat = fCMesh.FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat);
    
    pBC->Val1()(0,0) = -fFluidPressure;//-29.3;
    pBC->Val1()(1,1) = -fFluidPressure;//-29.3;

    
}


#include "pzgengrid.h"
#include "TPZVTKGeoMesh.h"

void TPZWellBoreAnalysis::TConfig::CreateMesh()
{
    if ((fGMesh.NElements() != 0) || (fGMesh.NNodes() != 0)) DebugStop();

    TPZManVector<REAL,3> x0(3,fInnerRadius),x1(3,fOuterRadius);
    x0[1] = 0;
    x1[1] = M_PI_2;
    x0[2] = 0.;
    x1[2] = 0.;
//    int numdiv = 20;
    TPZManVector<int,2> nx(fNx);
//    nx[0] = 50;
    TPZGenGrid gengrid(nx,x0,x1);
    REAL minsize = fDelx;//fCurrentConfig.fInnerRadius*M_PI_2/numdiv;
    REAL domainsize = fOuterRadius-fInnerRadius;
    REAL geoprogression = gengrid.GeometricProgression(minsize, domainsize, nx[0]);
    TPZManVector<REAL,2> geoprogressionvec(2,1.);
    geoprogressionvec[0] = geoprogression;
    gengrid.SetGeometricProgression(geoprogressionvec);
    gengrid.Read(&fGMesh);
    /// bottom side
    gengrid.SetBC(&fGMesh, 4, -4);
    /// outer side
    gengrid.SetBC(&fGMesh, 5, -3);
    /// left side
    gengrid.SetBC(&fGMesh, 6, -5);
    /// inner radius
    gengrid.SetBC(&fGMesh, 7, -2);
    
    /// wrap the mesh around
    int nnodes = fGMesh.NNodes();
    for (int in=0; in<nnodes; in++) {
        TPZManVector<REAL,3> cyl(3,0.), cart(3,0.);
        fGMesh.NodeVec()[in].GetCoordinates(cyl);
        cart[0] = cyl[0]*cos(cyl[1]);
        cart[1] = cyl[0]*sin(cyl[1]);
        fGMesh.NodeVec()[in].SetCoord(cart);
    }
#ifdef DEBUG
    std::ofstream out("Wellbore.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(&fGMesh, out,true);
#endif
    
    if (fGreater.size() == 0) {
        return;
    }
#ifdef LOG4CXX
    {
        std::stringstream sout;
        fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    // project the nodes on the elliptic boundaries
    for (int iy=0; iy < nx[1]+1; iy++) {
        int index = iy*(nx[0]+1);
        TPZManVector<REAL,3> co(3);
        fGMesh.NodeVec()[index].GetCoordinates(co);
        REAL xinit = co[0];
        bool changed = ProjectNode(co);
        if (changed) {
            REAL xadjust = co[0];
            int endindex = index+nx[0];
            TPZManVector<REAL,3> endco(3);
            fGMesh.NodeVec()[endindex].GetCoordinates(endco);
            REAL xend = endco[0];
            for (int ix=0; ix< nx[0]+1; ix++) {
                TPZManVector<REAL,3> nodeco(3);
                fGMesh.NodeVec()[index+ix].GetCoordinates(nodeco);
                REAL factor = (nodeco[0]-xinit)/(xend-xinit);
                REAL correctx = xadjust + factor*(xend-xadjust);
                nodeco[0] = correctx;
                fGMesh.NodeVec()[index+ix].SetCoord(nodeco);
            }
        }
    }
    // for the refined elements, put the nodes in the middle
    /*
    int nel = fGMesh.NElements();
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh.ElementVec()[el];
        if(!gel) continue;
        if(! gel->HasSubElement()) continue;
        for (int is = gel->NCornerNodes(); is < gel->NSides(); is++) {
            TPZManVector<int,5> midsidenodeindices;
            gel->MidSideNodeIndices(is, midsidenodeindices);
            if (midsidenodeindices.size() != 1) {
                DebugStop();
            }
            int nsnode = gel->NSideNodes(is);
            TPZManVector<REAL,3> xmid(3,0.);
            for (int isn=0; isn<nsnode; isn++) {
                TPZGeoNode *ptr = gel->SideNodePtr(is, isn);
                TPZManVector<REAL,3> xo(3);
                ptr->GetCoordinates(xo);
                for (int co=0; co<3; co++) xmid[co] += xo[co];
            }
            for(int co=0; co<3; co++) xmid[co] /= nsnode;
            TPZGeoNode *ptr = gel->NodePtr(midsidenodeindices[0]);
            ptr->SetCoord(xmid);
        }
    }
    for (int el=0; el<nel; el++) {
        TPZGeoEl *gel = fGMesh.ElementVec()[el];
        if (!gel || gel->MaterialId() != -2) continue;
        int nnodes = gel->NNodes();
        for (int in=0; in<nnodes; in++) {
            TPZGeoNode *ptr = gel->NodePtr(in);
            TPZManVector<REAL, 3> co(3);
            ptr->GetCoordinates(co);
            bool project = ProjectNode(co);
            if (project) {
                ptr->SetCoord(co);
            }
        }
    }
     */
}

void TPZWellBoreAnalysis::TConfig::ModifyWellElementsToQuadratic()
{
//    Phil, este eh o metodo que fiz!
//    Parece que a chamada ProjectNode(nodeCoord) (linha 1736) nao estah projetando legal!!!
//    Abraço.
//    
//    Caju.
    
    int nelem = fGMesh.NElements();
    for(int el = 0; el < nelem; el++)
    {
        TPZGeoEl * bcGel = fGMesh.ElementVec()[el];
        
#ifdef DEBUG
        if(!bcGel)
        {
            DebugStop();
        }
#endif
        
        if(bcGel->MaterialId() != -2)
        {
            continue;
        }
        
#ifdef DEBUG
        if(bcGel->Dimension() != 1)
        {
            DebugStop();
        }
#endif
        if(bcGel->HasSubElement() || bcGel->Father())
        {
            DebugStop();//este metodo nao funciona para elementos refinados e/ou filhos
        }
        int inner1DSide = 2;
        
        bcGel = TPZChangeEl::ChangeToQuadratic(&fGMesh, el);
        
        TPZGeoElSide bcSide(bcGel,inner1DSide);
        int nodeIdlocal = 2;
        int nodeIdnodevec = -1;
        nodeIdnodevec = bcGel->NodeIndex(nodeIdlocal);
        
        TPZManVector<REAL,3> nodeCoord(3,0.);
        fGMesh.NodeVec()[nodeIdnodevec].GetCoordinates(nodeCoord);
        ProjectNode(nodeCoord);
        
        TPZGeoElSide neighSide(bcSide.Neighbour());
        while(neighSide != bcSide)
        {
            TPZGeoEl * quadraticGel = TPZChangeEl::ChangeToQuadratic(&fGMesh, neighSide.Element()->Index());
            neighSide = quadraticGel->Neighbour(neighSide.Side());
        }
        
        fGMesh.NodeVec()[nodeIdnodevec].SetCoord(nodeCoord);
    }
#ifdef LOG4CXX
    if(logger->isDebugEnabled())
    {
        std::stringstream sout;
        fGMesh.Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
}

/// project the node on the boundary
// returns true if the coordinate was changed
bool TPZWellBoreAnalysis::TConfig::ProjectNode(TPZVec<REAL> &co)
{
    bool wasadjusted = false;
    for (int ellips = fGreater.size()-1; ellips >=0; ellips--) {
        if (co[1] > fSmaller[ellips]) {
            continue;
        }
        REAL a = fGreater[ellips];
        REAL b = fSmaller[ellips];
        REAL xadjust = a*sqrt(1.-co[1]*co[1]/(b*b));
        if (xadjust > co[0]) {
            co[0] = xadjust;
            wasadjusted = true;
        }
        // only one adjustment can be applied
    }
    return wasadjusted;
}

/// return the largest y-coordinate belonging to the last ellips
// this value will be used to identify the meaningful points to match the next ellips
REAL TPZWellBoreAnalysis::TConfig::MaxYfromLastBreakout()
{
    if(fGreater.size() == 0) return this->fInnerRadius;
    // we will use regula falsi
    int numell = fGreater.size();
    REAL lasta = fGreater[numell-1];
    REAL lastb = fSmaller[numell-1];
    REAL y = 0.;
    REAL dely = lastb/5;
    TPZManVector<REAL,2> co(2);
    while (fabs(dely) > lastb/1000.) {
        co[1] = y;
        co[0] = sqrt(fInnerRadius*fInnerRadius-y*y);
        REAL lastx = lasta*sqrt(1.-co[1]*co[1]/(lastb*lastb));
        ProjectNode(co);
        if (co[0] == lastx) {
            y += dely;
        }
        else
        {
            y -= dely;
            dely /= 10.;
        }
    }
    return y;
}



/// Initialize the Sandler DiMaggio object and create the computational mesh
void TPZWellBoreAnalysis::TConfig::CreateComputationalMesh(int porder)
{
    if ((fCMesh.NElements() != 0) || (fCMesh.NConnects() != 0)) DebugStop();

    fCMesh.SetReference(&fGMesh);
    int defaultporder = porder;
    fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &fCMesh;
    
#ifdef PV
    
    int materialid=1;
    bool planestrain=true;
    TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> &SD =fSDPV;
    TPZMatElastoPlastic2D<TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> > *PlasticSD = new TPZMatElastoPlastic2D< TPZPlasticStepPV<TPZSandlerExtended,TPZElasticResponse> >(materialid,planestrain);
    
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    TPZTensor<REAL> initstress(0.),finalstress(0.);
    REAL hydro = fConfinement.I1();
    hydro -= SD.fYC.fA*SD.fYC.fR;
    hydro /= 3.;
    finalstress.XX() = hydro;
    finalstress.YY() = hydro;
    finalstress.ZZ() = hydro;
    
    PrepareInitialMat(SD, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = fConfinement;
    PrepareInitialMat(SD, initstress, finalstress, 10);
    SD.ResetPlasticMem();
    
    TPZMaterial *plastic(PlasticSD);
    
    PlasticSD->SetPlasticity(SD);
    compmesh1->InsertMaterialObject(plastic);
    
    CmeshWell(compmesh1, plastic, fConfinement, fFluidPressure);
    compmesh1->AutoBuild();
#else
    //    //TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2>::UncDeepSandTest(obj.fCurrentConfig.fSD);
    TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> &SD = fSD;
    
    SD.SetResidualTolerance(1.e-10);
    SD.fIntegrTol = 10.;
    int materialid = 1;
    bool planestrain = true;
    TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> > *PlasticSD = new TPZMatElastoPlastic2D<TPZSandlerDimaggio<SANDLERDIMAGGIOSTEP2> >(materialid,planestrain);
    
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(compmesh1);
    
    
    
    
    TPZTensor<REAL> initstress,finalstress;
    REAL hydro = fConfinement.I1();
    hydro -= SD.fYC.fA*SD.fYC.fR;
    hydro /= 3.;
    finalstress.XX() = hydro;
    finalstress.YY() = hydro;
    finalstress.ZZ() = hydro;
    
    PrepareInitialMat(SD, initstress, finalstress, 10);
    initstress = finalstress;
    finalstress = fConfinement;
    PrepareInitialMat(SD, initstress, finalstress, 10);
    SD.ResetPlasticMem();
    
    TPZMaterial *plastic(PlasticSD);
    
    PlasticSD->SetPlasticity(SD);
    compmesh1->InsertMaterialObject(plastic);
    
    CmeshWell(compmesh1, plastic, fConfinement, fFluidPressure);
    compmesh1->AutoBuild();
    
#endif
    
    
    
}

#include "pzelasmat.h"

/// Configure the wellbore analysis to perform a linear analysis
void TPZWellBoreAnalysis::LinearConfiguration(int porder)
{
    fCurrentConfig.fCMesh.SetReference(&fCurrentConfig.fGMesh);
    int defaultporder = porder;
    fCurrentConfig.fCMesh.SetDefaultOrder(defaultporder);
    TPZCompMesh *compmesh1 = &fCurrentConfig.fCMesh;
    TPZElasticityMaterial *elasmat = new TPZElasticityMaterial(1);
    ConfigureLinearMaterial(*elasmat);
    compmesh1->InsertMaterialObject(elasmat);
    CmeshWell(compmesh1, elasmat, fCurrentConfig.fConfinement, fCurrentConfig.fFluidPressure);
    compmesh1->AutoBuild();
    
}


/// Set the parameters of the linear material
void TPZWellBoreAnalysis::ConfigureLinearMaterial(TPZElasticityMaterial &mat)
{
    REAL E,nu, lambda,G;
    G = fCurrentConfig.fSD.fER.fMu;
    lambda = fCurrentConfig.fSD.fER.fLambda;
    E=G*(3.*lambda+2.*G)/(lambda+G);
    nu = lambda/(2.*(lambda+G));
    mat.SetElasticity(E, nu);
    
    REAL sigxx = fCurrentConfig.fConfinement.XX();
    REAL sigyy = fCurrentConfig.fConfinement.YY();
    REAL sigxy = fCurrentConfig.fConfinement.XY();
    REAL sigzz = fCurrentConfig.fConfinement.ZZ();
    mat.SetPreStress(sigxx, sigyy, sigxy, sigzz);
    
}

/// Compute the linear elastic stiffness matrix
void TPZWellBoreAnalysis::ComputeLinearMatrix()
{
    TPZSkylineStructMatrix skylstr(&fCurrentConfig.fCMesh);
    skylstr.SetNumThreads(8);
    TPZFMatrix<STATE> rhs;
    TPZCompMesh *compmesh1 = &fCurrentConfig.fCMesh;

    TPZMaterial *keepmat = compmesh1->MaterialVec()[1];
    TPZElasticityMaterial *elasmat = new TPZElasticityMaterial(1);
    ConfigureLinearMaterial(*elasmat);
    compmesh1->InsertMaterialObject(elasmat);
    fLinearMatrix = skylstr.CreateAssemble(rhs, NULL);
    compmesh1->InsertMaterialObject(keepmat);
    delete elasmat;
}

STATE TPZWellBoreAnalysis::TConfig::ComputeFarFieldWork()
{
    STATE work=0.,nx=0.,ny=0.,radius=0.;
    TPZManVector<REAL,3> ksi(3,0.),sol(3,0.),xvec(3,0.);
    TPZFMatrix<STATE> jac,jacinv,axes,invjac;
    STATE detjac;
    
    int nel = fCMesh.NElements();
    fCMesh.Reference()->ResetReference();
    fPostprocess.Mesh()->LoadReferences();
    
    for (int el=0; el<nel; el++)
    {
        
        TPZCompEl * cel = fCMesh.ElementVec()[el];
        if(!cel)
        {
            continue;
        }
        TPZGeoEl *gel = cel->Reference();
        if(!gel)
        {
            DebugStop();
        }
        int matid = gel->MaterialId();
        if(matid!=-3)
        {
            continue;
        }
        
        TPZGeoElSide gelside(gel,2);
        TPZStack<TPZCompElSide> stackcompelside;
        gelside.EqualLevelCompElementList(stackcompelside, 0, 0);
        if(stackcompelside.size()!=1)
        {
            DebugStop();
        }
        
        TPZTransform t1(1);
        TPZCompElSide compneigh = stackcompelside[0];
        TPZGeoElSide neighbour = stackcompelside[0].Reference();
        gelside.SideTransform3(neighbour,t1);
        int sidefrom,sideto = neighbour.Element()->NSides()-1;
        sidefrom=neighbour.Side();
        TPZTransform t2 = neighbour.Element()->SideToSideTransform(sidefrom, sideto);
        
        
        TPZInterpolationSpace *intel = dynamic_cast<TPZInterpolationSpace *>(cel);
        const TPZIntPoints &rule = intel->GetIntegrationRule();
        REAL weight;
        int npoints = rule.NPoints();
        
        for(int i=0;i<npoints;i++)
        {
            rule.Point(i, ksi, weight);
            
            TPZManVector<REAL,3> ksi1(1,0.), ksi2(2,0.);
            t1.Apply(ksi, ksi1);
            t2.Apply(ksi1, ksi2);
            compneigh.Element()->Solution(ksi2, 30, sol);
            
            gel->X(ksi, xvec);
            intel->GetIntegrationRule().Point(0, ksi, weight);
            radius=sqrt(xvec[0]*xvec[0]+xvec[1]*xvec[1]);
            nx=xvec[0]/radius;
            ny=xvec[1]/radius;
            gel->Jacobian(xvec,jac,axes,detjac,invjac);
            work +=(fConfinement.XX()*nx*sol[0]+fConfinement.YY()*ny*sol[1])*weight*fabs(detjac);
            
        }
        
    }
    cout << "\n WORK "<< work << endl;
    return work;
    
}

void TPZWellBoreAnalysis::Print(std::ostream &out)
{
    out << "Number of configurations " << fSequence.size() << std::endl;
    out << "fCurrentConfig \n";
    fCurrentConfig.Print(out);
    out << "Other configurations \n";
    std::list<TPZWellBoreAnalysis::TConfig>::iterator it;
    for (it= fSequence.begin(); it != fSequence.end(); it++) {
        (*it).Print(out);
    }
    if (fLinearMatrix) {
        fLinearMatrix->Print("Linear Matrix",out);
    }
    else
    {
        out << "No linear matrix ";
    }
}

void TPZWellBoreAnalysis::TConfig::Print(ostream &out)
{
    out << "fInnerRadius " << fInnerRadius << endl;
    out << "fOuterRadius " << fOuterRadius << endl;
    out << "fNx " << fNx << endl;
    out << "fDelx " << fDelx << endl;
    out << "fGreater " << fGreater << endl;
    out << "fSmaller " << fSmaller << endl;
    out << "fConfinement " << fConfinement << endl;
    out << "fSD ";
    fSD.Print(out);
    out << "fFluidPressure " << fFluidPressure << endl;
    out << "fGMesh ";
    fGMesh.Print(out);
    out << "fCMesh ";
    fCMesh.Print(out);
    fAllSol.Print("fAllSol = ",out);
    out << "fPlasticDeform " << fPlasticDeformSqJ2 << endl;
    fPostprocess.Print("fPostProcess", out);
}
