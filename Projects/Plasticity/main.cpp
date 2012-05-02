//$Id: poroelastoplastic.cpp,v 1.53 2010-06-11 22:13:02 diogo Exp $


/***************************************************************************
 *   Copyright (C) 2008 by Erick Slis   *
 *      *
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 *   This program is distributed in the hope that it will be useful,       *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         *
 *   GNU General Public License for more details.                          *
 *                                                                         *
 *   You should have received a copy of the GNU General Public License     *
 *   along with this program; if not, write to the                         *
 *   Free Software Foundation, Inc.,                                       *
 *   59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.             *
 ***************************************************************************/


#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
/*
 using namespace std;
 
 int main(int argc, char *argv[])
 {
 cout << "Hello, world!" << endl;
 
 return EXIT_SUCCESS;
 }
 */
#include "pzelasmat.h"
#include "BrazilianTestGeoMesh.h"
#include "pzelastoplastic2D.h"
#include "pzelctemp.h" // TPZIntelGen<TSHAPE>
#include "pzshapecube.h" // TPZShapeCube
#include "pzcompelwithmem.h"
#include "pzelastoplastic.h"
#include "pzporous.h"
#include "TPZLadeKim.h"
#include "TPZSandlerDimaggio.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzelastoplasticanalysis.h"
#include "pzporoanalysis.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"

#include "TPZTensor.h"

#include "pzelast3d.h"
#include "pzbndcond.h"

#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"

#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include <sstream>

using namespace pzshape; // needed for TPZShapeCube and related classes

#include "pzlog.h"

#ifdef LOG4CXX // LOG4CXX may be defined alone or with LOG4CXX_PLASTICITY. The latter shall not be used alone.
#include <log4cxx/logger.h>
#include <log4cxx/basicconfigurator.h>
#include <log4cxx/propertyconfigurator.h>
#endif

#ifdef LOG4CXX
static LoggerPtr PEPMainlogger(Logger::getLogger("poroelastoplastic.main"));
#endif


void InitializeLOG()
{
#ifdef LOG4CXX
    
    std::string path;
    std::string configfile;
#ifdef HAVE_CONFIG_H
    path = POROELASTOPLASTICSOURCEDIR;
    path += "/src/";
    cout << path.c_str() << endl;
    cout.flush();
#else
    path = "";
#endif
    configfile = path;
    
    configfile += "log4cxx.cfg";
    log4cxx::PropertyConfigurator::configure(configfile.c_str());
    
    std::stringstream sout;
    sout << __PRETTY_FUNCTION__ << "\nLOG4CXX configured.\n"
    << "LOG4CXX config file:" << configfile;
    LOGPZ_INFO(PEPMainlogger,sout.str().c_str());
    
#endif
}

void SolveSistII(TPZAnalysis &an, TPZCompMesh *fCmesh);
void CMeshGid(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);

void CMeshTwoMaterials(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> matelast,TPZAutoPointer<TPZMaterial> matplastic);

void CMesh2DBCII(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat);
void SetUPPostProcessVariablesII(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );
void ManageIterativeProcessII(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
							int BCId,int BCId2, int nsteps, REAL PGRatio,
							TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
							TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
							TPZPostProcAnalysis * ppAnalysis, int res);

void ManageIterativeProcessIII(std::ostream &out,REAL tol,int numiter,
                               int BCId, int nsteps, REAL PGRatio,
                               TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
                               TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
                               TPZPostProcAnalysis * ppAnalysis, int res);


void BrazilianPlasticAnalysis2D();
void SetUPPostProcessElasticVariables2D(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames );


void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
                            int BCId,int BCId2, int nsteps, REAL PGRatio,
                            TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
                            TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End)
{
    
    if(!analysis.Mesh())return;
    
    // computing the initial value for the PG progression such that its sum equals one;
    REAL a0;
    if(fabs(PGRatio - 1.) < 1.e-3)
    {
        a0 = 1. / REAL(nsteps);
    }else{
        a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
    }
    TPZFNMatrix<36,REAL > val1(6,6,0.), deltaVal1(6,6,0.);
    TPZFNMatrix< 6,REAL > val2(6,1,0.), deltaVal2(6,1,0.);
    
    deltaVal1 = val1End;
    deltaVal1.ZAXPY(-1., val1Begin);
    deltaVal2 = val2End;
    deltaVal2.ZAXPY(-1., val2Begin);
    
    // ZAXPY operation: *this += alpha * p			
    
    TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
    TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
    if(!pBC)return;
    
    int i;
    for(i = 0; i < nsteps; i++)
    {
        REAL stepLen;
        if(fabs(PGRatio - 1.) < 1.e-3)
        {
            stepLen = REAL(i+1) / REAL(nsteps);
        }else{
            stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
        }
        
        val1 = val1Begin;
        val1.ZAXPY(stepLen, deltaVal1);
        val2 = val2Begin;
        val2.ZAXPY(stepLen, deltaVal2);
        
        pBC->Val1() = val1;
        pBC->Val2() = val2;
        
        analysis.IterativeProcess(out, tol, numiter);
        
        
        //analysis.AcceptSolution();
        
    }
    
}



void SolveSistLin2(TPZAnalysis &an, TPZCompMesh *fCmesh)
{
    //TPZSkylineStructMatrix sky(fCmesh);
    //	an.SetStructuralMatrix(sky);
    //	TPZStepSolver step;
    //	step.SetDirect(ECholesky);
    //	an.SetSolver(step);
    //	an.SetCompMesh(fCmesh);
    //	an.Run(cout);
    
    
    TPZSkylineStructMatrix full(fCmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver<REAL> step;
    //full.SetNumThreads(8);
    step.SetDirect(ELDLt);
    an.SetSolver(step);
}

template <class T>
void WellboreLoadTestWithoutPostProcess(stringstream & fileName, T & mat, 
                                        REAL loadMultipl, REAL plasticTol)
{
    REAL L, LDA, alpha, GammaSea, MudWeight, kh, kH, OCR, s;
    s = loadMultipl;
    REAL SigmaV=0., SigmaH=0., Sigmah=0.;
    REAL Pa = 14.7;
    int ncirc, ioRatio, pOrder, valType;
    
    cout << "\nMesh data: ncirc?(int) ";
    //ncirc = 10;
    ncirc=10;//cin >> ncirc;
    
    cout << "Mesh data: ioratio?(sugg. 10.) ";
    //ioRatio =  10.;
    ioRatio=10.;//cin >> ioRatio;
    
    cout << "Mesh data: pOrder? ";
    //pOrder =2;
    pOrder=3;//cin >> pOrder;
    
    fileName << "_ncirc" << ncirc << "_IO" << ioRatio << "_p" << pOrder;
    
    cout << "\nSelect one of the following load case (Load):";
    cout << "\n0) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.8 kH=0.9 OCR=1.10";
    cout << "\n1) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.6 kH=0.8 OCR=1.10";
    cout << "\n2) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=0.5";
    cout << "\n3) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=1.0";
    cout << "\n4) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=9.5 kh=0.80 kH=1.30 OCR=1.10 Biot=1.0";
    cout << "\n";
    
    //valType = 0 ;
    valType=0;//cin >> valType;
    
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
    
    
    TPZCompMesh * pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
    
    //building analysis
    TPZElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
    
    SolveSistLin2(EPAnalysis,pCMesh);
    EPAnalysis.SetBiCGStab(5000, 1.e-12);
    
    //    void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
    //                                int BCId,int BCId2, int nsteps, REAL PGRatio,
    //                                TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
    //                                TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End)
    REAL tol = 1.e-5;
    int numiter = 30;
    int nsteps = 1;
    REAL pgratio = 1.;
    
    
    
    ManageIterativeProcess(EPAnalysis,cout,tol, numiter,
                           -7 /*BCId*/, 2 /*nsteps*/, nsteps, pgratio,
                           BeginStress/*val1Begin*/, EndStress/*val1End*/,
                           val1/*val2Begin*/, val1/*val2End*/);
    
    ManageIterativeProcess(EPAnalysis,cout, tol, numiter,
                           -7 /*BCId*/, 2 /*nsteps*/, nsteps,pgratio,
                           EndStress/*val1Begin*/, EndStress2/*val1End*/,
                           val1/*val2Begin*/, val1/*val2End*/);
    
    cout << "\n Problem solved ";
    
    
    return;
    
    
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
    //ncirc = 10;
    ncirc=8;//cin >> ncirc;
    
    cout << "Mesh data: ioratio?(sugg. 10.) ";
    //ioRatio =  10.;
    ioRatio=10.;//cin >> ioRatio;
    
    cout << "Mesh data: pOrder? ";
    //pOrder =2;
    pOrder=2;//cin >> pOrder;
    
    fileName << "_ncirc" << ncirc << "_IO" << ioRatio << "_p" << pOrder;
    
    cout << "\nSelect one of the following load case (Load):";
    cout << "\n0) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.8 kH=0.9 OCR=1.10";
    cout << "\n1) L=3000 LDA=1200 alpha=0.8 GammaSea=8.6 MudWeight=9.2 kh=0.6 kH=0.8 OCR=1.10";
    cout << "\n2) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=0.5";
    cout << "\n3) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=10. kh=0.96 kH=1.09 OCR=1.10 Biot=1.0";
    cout << "\n4) L=5227 LDA=2135 alpha=0.5 GammaSea=9.5 MudWeight=9.5 kh=0.80 kH=1.30 OCR=1.10 Biot=1.0";
    cout << "\n";
    
    //valType = 0 ;
    valType=0;//cin >> valType;
    
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
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Material overconsolidated at " << OCStress;
        sout << "\nAnd reset to the initial value of " << FarFieldStress;
        sout << "\nInitial Material settings:\n:";
        EPMat.Print(sout, 1/*Memory*/);
        LOGPZ_INFO(PEPMainlogger,sout.str().c_str());
    }
#endif
    
    //End of material initialization
    
    
    TPZCompMesh * pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
    
    //TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
    
    //building analysis
    TPZElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
    
    SolveSistLin2(EPAnalysis,pCMesh);
    EPAnalysis.SetBiCGStab(5000, 1.e-12);
    
    // Preparing Post Process
    
    TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
    
    TPZSkylineStructMatrix structmatrix(PPAnalysis.Mesh());
    PPAnalysis.SetStructuralMatrix(structmatrix);
    
    EPAnalysis.TransferSolution(PPAnalysis);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZVec<std::string> PostProcVars, scalNames, vecNames;
    
    scalNames.Resize(7);
    scalNames[0] = "Alpha";
    scalNames[1] = "PlasticSteps";
    scalNames[2] = "VolElasticStrain";
    scalNames[3] = "VolPlasticStrain";
    scalNames[4] = "VolTotalStrain";
    scalNames[5] = "I1Stress";
    scalNames[6] = "J2Stress";
    
    vecNames.Resize(5);
    vecNames[0] = "Displacement";
    vecNames[1] = "NormalStress";
    vecNames[2] = "ShearStress";
    vecNames[3] = "NormalStrain";
    vecNames[4] = "ShearStrain";
    
    PostProcVars.Resize(scalNames.NElements()+vecNames.NElements());
    int i, k=0;
    for(i = 0; i < scalNames.NElements(); i++)
    {
        PostProcVars[k] = scalNames[i];
        k++;
    }
    for(i = 0; i < vecNames.NElements(); i++)
    {
        PostProcVars[k] = vecNames[i];
        k++;
    }
    
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
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "\nend Material settings:\n:";
        EPMat.Print(sout, 1/*Memory*/);
        LOGPZ_INFO(PEPMainlogger,sout.str().c_str());
    }
#endif
    
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
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "Material overconsolidated at " << OCStress;
        sout << "\nAnd reset to the initial value of " << FarFieldStress;
        sout << "\nInitial Material settings:\n:";
        EPMat.Print(sout, 1/*Memory*/);
        LOGPZ_INFO(PEPMainlogger,sout.str().c_str());
    }
#endif
    
    //End of material initialization
    
    
    TPZCompMesh * pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
    
    TPZPoroElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh); // self explanatory
    
    {// Replacing BC to include Pore Pressure
        TPZAutoPointer<TPZMaterial> mat(&EPMat);
        TPZFMatrix<REAL> val2(4,1,0.), val1(4,4,0.);
        val2.Zero();
        val2(0,0) = 1.;
        val2(1,0) = 1.;
        val2(3,0) = 1.;
        val1(3,3) = 0;// no pore pressure change
        TPZAutoPointer<TPZMaterial> bc = mat->CreateBC(mat, -5, 3, val1, val2);//Directional 
        pCMesh->InsertMaterialObject(bc);
    }
    
    //building analysis
    TPZPoroElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
    
    EPAnalysis.SetBiCGStab(50000, 1.e-10);
    
    // Preparing Post Process
    
    TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
    
    TPZVec<int> PostProcMatIds(1,1);
    TPZVec<std::string> PostProcVars, scalNames, vecNames;
    
    scalNames.Resize(8);
    scalNames[0] = "PorePressure";
    scalNames[1] = "Alpha";
    scalNames[2] = "PlasticSteps";
    scalNames[3] = "VolElasticStrain";
    scalNames[4] = "VolPlasticStrain";
    scalNames[5] = "VolTotalStrain";
    scalNames[6] = "I1Stress";
    scalNames[7] = "J2Stress";
    
    vecNames.Resize(5);
    vecNames[0] = "Displacement";
    vecNames[1] = "NormalStress";
    vecNames[2] = "ShearStress";
    vecNames[3] = "NormalStrain";
    vecNames[4] = "ShearStrain";
    
    PostProcVars.Resize(scalNames.NElements()+vecNames.NElements());
    int i, k=0;
    for(i = 0; i < scalNames.NElements(); i++)
    {
        PostProcVars[k] = scalNames[i];
        k++;
    }
    for(i = 0; i < vecNames.NElements(); i++)
    {
        PostProcVars[k] = vecNames[i];
        k++;
    }
    
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
    for(i = 0; i < 10; i++)
    {
        REAL currentTime = pow(10.,i);
        deltaT = currentTime - time;
        time = currentTime;
        cout << "\n Evoluting to time " << currentTime;
        
        EPAnalysis.SetDeltaT(deltaT);
        
        EPAnalysis.Run(cout, 1.e-4, 10, &PPAnalysis, pOrder);
    }
    
    PPAnalysis.PostProcess(pOrder);
    
    cout << "\nProblem Solved. Accepting new Solutions\n";
    
    cout << "\nClosing Mesh\n";
    
    PPAnalysis.CloseGraphMesh();
    
    cout << "\nExiting\n";
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "\nend Material settings:\n:";
        EPMat.Print(sout, 1/*Memory*/);
        LOGPZ_INFO(PEPMainlogger,sout.str().c_str());
    }
#endif
    
    return;
    
    
}

#include "TPZPlasticityTest.h"
//
void MaterialPointTests()
{
    
    cout << "\nChoose Plasticity test:";
    cout << "\n0 - Isotropic compression ";
    cout << "\n1 - Biaxial Tests ";
    cout << "\n2 - Uniaxial traction ";
    cout << "\n";
    int choice = 1;
   //cin >> choice;
    
    switch(choice)
    {
        case(0):
            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
            cout << "\n0 - Lade - Kim ";
            cout << "\n1 - Sandler Dimaggio ";
            cout << "\n2 - Drucker Prager ";
            cout << "\n";
            int choice2;
           // cin >> choice2;
            choice2 = 0;
            switch(choice2)
        {
            case(0):
                LKIsotropicCompression();
                break;
            case(1):
                SandlerDimaggioIsotropicCompression();
                break;
            case(2):
                DruckerIsotropicCompression();
                break;
        }
            
            
            break;
        case(1):
            LKBiaxialTest();
            break;
        case(2):
            cout << "NOT IMPLEMENTED YET";
            //    LadeKim_ReversalTest();
            break;
        default:
            cout << "Unknown Test Type. Exiting...";
    }
    
}


//#include "BrazilianTestGeoMesh.h"

#include "fad.h"

Fad<double> func(const Fad<double>& x, const Fad<double>& y)
{
//	double z=sqrt(x);
    /*
     //Mathematica
    f = {y Sqrt[x] + Sqrt[x]}
    f /. {x -> 1, y -> 2}
    dfx = D[f, x];
    dfy = D[f, y];
    dfx /. {x -> 1., y -> 2}
    dfy /. {x -> 1., y -> 2}
     */
	return y*sqrt(x)+sqrt(x);
}



int main()
{
    
    
    InitializeLOG();
    TPZLadeKim LK;
    LK.PlainConcrete(LK);
    TPZDruckerPrager DP;
   DP.PlainConcreteMPa(DP);
//    REAL dirMult = 60.;
//    MultiDirectionsMaterialPointTest(LK,dirMult);
  
    
//    double x,y,f;     // Declare variables x,y,f
//    x=1;                 // Initialize variable x
//    x.diff(0,2);         // Differentiate with respect to x (index 0 of 2)
//    y=2;                 // Initialize variable y
//    y.diff(1,2);         // Differentiate with respect to y (index 1 of 2)
//    f=func(x,y);         // Evaluate function and derivatives
//    double fval=f.x();   // Value of function
//    double dfdx=f.d(0);  // Value of df/dx (index 0 of 2)
//    double dfdy=f.d(1);  // Value of df/dy (index 1 of 2)
    
//    cout << "f(x,y)=" << fval << endl;
//    cout << "df/dx(x,y)=" << dfdx << endl;
//    cout << "df/dy(x,y)=" << dfdy << endl;

    
    Fad<REAL> x,y,f;
    x=1.;
    x.diff(0,2);
    y=2.;
    y.diff(1,2);
    f = func(x,y);
    double fval=f.val();   // Value of function
    double dfdx=f.d(0);  // Value of df/dx (index 0 of 2)
    double dfdy=f.d(1);  // Value of df/dy (index 1 of 2)
    
    cout << "f(x,y)=" << fval << endl;
    cout << "df/dx(x,y)=" << dfdx << endl;
    cout << "df/dy(x,y)=" << dfdy << endl;
    

    
    
//  TPZPlasticTest::DruckerTest();
 //   LadeKimTriaxialLooseSand();
  // LKIsotropicCompression();
  //  LKBiaxialTest();
  //  SandlerDimaggioIsotropicCompression();
   // TPZTensor<REAL> sigma,eps;
   //sigma.XX()= 435.286;
   // sigma.YY()=170.146;
   // sigma.ZZ()= -327.067;
    //TPZLadeKim LK;
   // LK.PlainConcrete(LK);
   // LK.ApplyLoad(sigma,eps);
  //  FragGranade();
   // LKIsotropicCompression();

  //  BrazilianPlasticAnalysis2D();
    
  //  TPZGeoMesh * mesh = new TPZGeoMesh;
  //  mesh = BrazilianTestGeoMesh::MalhaPredio();
//    mesh = MalhaPredio();
/*    

    MultiDirectionsMaterialPointTest(DP,dirMult);
    cout << "\nRuning finished " << endl;
 */
  
  //  TPZPlasticTest::LoadTest("testeFineSilicaIsotropic.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaIsotropic0.001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaIsotropic0.0001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaIsotropic0.00001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaIsotropic0.000001.loadpath");
 //   TPZPlasticTest::LoadTest("testeFineSilicaOCR1_0.0001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaOCR1_0.00001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaOCR1.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaOCR10_0.0001.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaOCR10_0.00001.loadpath");
   // TPZPlasticTest::LoadTest("testeFineSilicaOCR10.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaPEP.loadpath");
  //  TPZPlasticTest::LoadTest("testeFineSilicaPEP2.loadpath");
 //   TPZPlasticTest::LoadTest("testeLooseSacrRiverSandHidrostatico0.001.loadpath");
//    TPZPlasticTest::LoadTest("testeMcCormicRanchSand_HID.loadpath");
//    TPZPlasticTest::LoadTest("testeMcCormicRanchSand_HIDstrain.loadpath");
//    TPZPlasticTest::LoadTest("testeMcCormicRanchSand_PLT_0.6_0.0001.loadpath");
   // TPZPlasticTest::LoadTest("testeMcCormicRanchSand_PLT_0.6_0.00001.loadpath");
 
    
    TPZPlasticTest::PlasticIntegratorCheck(LK);
  //  TPZPlasticTest::LoadTest("testeLooseSacrRiverSandHidrostatico0.001.loadpath");
       cout << "\nRuning finished " << endl;
//    TPZPlasticTest::LoadTest("testeMcCormicRanchSand_PLT_0.6_0.00001.loadpath");
//    TPZPlasticTest::LoadTest("testePlainStrain9b_0.00001.loadpath");
 //   TPZPlasticTest::LoadTest("DPPlainConcreteXX.loadpath");
  
  //  TPZPlasticTest::LoadTest("testeMcCormicRanchSand_HIDstrain.loadpath");
 
  //  LKIsotropicCompression();
  //  SandlerDimaggioIsotropicCompression();
  //   DruckerIsotropicCompression();
   // LKBiaxialTest();
   
/* 
    cout << "\nPlease enter test type:";
    cout << "\n0 Material Point test ";
    cout << "\n1 Finite Elements test ";
    
    int testnum;
    testnum =1;
   //100 cin >> testnum;
    
    switch (testnum) {
        case 0:
        {
            MaterialPointTests();
            break;
        }
        case 1:
        {
            int testNumber, matNumber;
            TPZPlasticBase *pMat;
            TPZLadeKim * pLK = NULL;
            TPZSandlerDimaggio * pSD = NULL;
            typedef TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> TPZDruckerPrager;
            TPZDruckerPrager * pDP = NULL;
            REAL loadMultipl;
            stringstream fileName, copyStr;
            REAL plasticTol;
            
            cout << "\nPlease enter test type:";
            cout << "\n0) Wellbore Drilling Load";
            cout << "\n1) Wellbore Drilling Load - Porous Medium";
            cout << "\n2) Wellbore Drilling Load # with out postprocess # (FOR TESTS PURPOSE) ";
            
            testNumber=2;//cin >> testNumber;
            
            cout << "\nMaterial Type:";
            cout << "\n0)Lade Kim: FineSilicaSand";
            cout << "\n1)Lade Kim: LooseSacrRiverSand";
            cout << "\n2)Lade Kim: DenseSacrRiverSand";
            cout << "\n3)Lade Kim: PlainConcrete";
            cout << "\n4)Sandler Dimaggio: McCormic Ranch Sand";
            cout << "\n5)Sandler Dimaggio: McCormic Ranch Sand Mod";
            cout << "\n6)Sandler Dimaggio: Unconsolidated Deep Reservoir Sandstone [psi]";
            cout << "\n7)Sandler Dimaggio: Unconsolidated Deep Reservoir Sandstone [MPa]";
            cout << "\n8)Sandler Dimaggio: PRSMat [MPa]";
            cout << "\n9)Drucker Prager (Inscr MC): PRSMat [MPa]";
            cout << "\n10)Drucker Prager (Circunscr MC): PRSMat [MPa]";
            cout << "\n";
            
            
            matNumber=8;//cin >> matNumber;
            
            switch(matNumber)
            {
                case(0):
                    pLK = new TPZLadeKim();
                    TPZLadeKim::FineSilicaSand(*pLK);
                    pMat = pLK;
                    fileName << "_LKFS";
                    loadMultipl = -1.;
                    break;
                case(1):
                    pLK = new TPZLadeKim();
                    TPZLadeKim::LooseSacrRiverSand(*pLK);
                    pMat = pLK;
                    fileName << "_LKLS";
                    loadMultipl = -1.;
                    break;
                case(2):
                    pLK = new TPZLadeKim();
                    TPZLadeKim::DenseSacrRiverSand(*pLK);
                    pMat = pLK;
                    fileName << "_LKDS";
                    loadMultipl = -1.;
                    break;
                case(3):
                    pLK = new TPZLadeKim();
                    TPZLadeKim::PlainConcrete(*pLK);
                    pMat = pLK;
                    fileName << "_LKPC";
                    loadMultipl = -1.;
                    break;
                case(4):
                    pSD = new TPZSandlerDimaggio();
                    TPZSandlerDimaggio::McCormicRanchSand(*pSD);
                    pMat = pSD;
                    fileName << "_SDMc";
                    loadMultipl = -0.001;
                    break;
                case(5):
                    pSD = new TPZSandlerDimaggio();
                    TPZSandlerDimaggio::McCormicRanchSandMod(*pSD);
                    pMat = pSD;
                    fileName << "_SDMM";
                    loadMultipl = -0.001;
                    break;
                case(6):
                    pSD = new TPZSandlerDimaggio();
                    TPZSandlerDimaggio::UncDeepSandResPSI(*pSD);
                    pMat = pSD;
                    fileName << "_SDDS";
                    loadMultipl = -1;
                    break;
                case(7):
                    pSD = new TPZSandlerDimaggio();
                    TPZSandlerDimaggio::UncDeepSandResMPa(*pSD);
                    pMat = pSD;
                    fileName << "_SDDSMPa";
                    loadMultipl = -1/145.03773801;
                    break;
                case(8):
                    pSD = new TPZSandlerDimaggio();
                    TPZSandlerDimaggio::PRSMatMPa(*pSD);
                    pMat = pSD;
                    fileName << "_PRSLMPa";
                    loadMultipl = -1/145.03773801;
                    break;
                case(9):
                    pDP = new TPZDruckerPrager();
                    pDP->fYC.SetUp( 29.7/180. *M_PI ,1);
                    pDP->fTFA.SetUp( 12.8, 1.);
                    pDP->fER.SetUp( 29269., 0.203);
                    pMat = pDP;
                    fileName << "_PRDPInscMPa";
                    loadMultipl = -1/145.03773801;
                    break;
                case(10):
                    pDP = new TPZDruckerPrager();
                    pDP->fYC.SetUp( 29.7/180. * M_PI ,0);
                    pDP->fTFA.SetUp( 12.8,1.);
                    pDP->fER.SetUp( 29269.,0.203);
                    pMat = pDP;
                    fileName << "_PRDPCircMPa";
                    loadMultipl = -1/145.03773801;
                    break;
                default:
                    cout << "\nUnhandled Test Type. Exiting...";
                    break;
            }
            
            cout << "\nPlastic Integration Tolerance:(sugg. 0.0001) ";
            
            plasticTol=0.0001;//cin >> plasticTol;
            
            fileName << "_pTol" << plasticTol;
            
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
                case(2):
                    copyStr << fileName.str();
                    fileName.str(""); // clearing the buffer
                    fileName << "PorousWB" << copyStr.str();
                    if(pLK)WellboreLoadTestWithoutPostProcess(fileName, *pLK, loadMultipl, plasticTol);
                    if(pSD)WellboreLoadTestWithoutPostProcess(fileName, *pSD, loadMultipl, plasticTol);
                    if(pDP)WellboreLoadTestWithoutPostProcess(fileName, *pDP, loadMultipl, plasticTol);
                    break;
                default:
                    cout << "\nUnhandled Test Type. Exiting...";
                    delete pMat;
            }
            break;
        }
        default:
        {
            cout << "\nUnhandled Material Type. Exiting...";
            break;
        }
    }
    
*/
    
    return EXIT_SUCCESS;
    
    
}

/*
void BrazilianElasticAnalysis()
{
	int h=4;
	int order = 2;
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	BrazilianTestGeoMesh::TransformBlendToLinearMesh(MESH,h);
	ofstream arg1("ElasticMaterialGeoMesh.txt");
	MESH->Print(arg1);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
	TPZVec<REAL> force(3,1);
	force[0]=0.;
	force[1]=0.;
	force[2]=0.;
	
	TPZElasticity3D *elast = new TPZElasticity3D(1,20000.,0.,force);
	
	CMESH->SetDefaultOrder(order);
	TPZAutoPointer<TPZMaterial> elastt = elast;
	CMESH->InsertMaterialObject(elastt);
	
	CMeshBC(CMESH,elastt);
	
	ofstream arg("CMESHELASTIC.txt");
	CMESH->Print(arg);
	
	TPZAnalysis an;
	TPZVec<int> PostProcMatIds(1,1);
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	
	SetUPPostProcessElasticVariables(PostProcVars,scalNames,vecNames);
	
	SolveSistLin(an,CMESH);
	
	cout << " \n ELASTIC-SOL \n";
	cout << an.Solution() << endl;
	an.DefineGraphMesh(3, scalNames,vecNames,"ELASTIC1.vtk");
	an.PostProcess(0);
	an.CloseGraphMesh();
	
}
*/

void SetUPPostProcessElasticVariables2D(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	vecnames.Resize(1);
	vecnames[0] = "displacement";
	
	scalnames.Resize(4);
	scalnames[0]="SigmaX";
	scalnames[1]="SigmaY";
	scalnames[2]="tau_xy";
	scalnames[3]="Pressure";
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


void PrepareInitialMat(TPZPlasticBase & mat, TPZTensor<REAL> &initialStress, TPZTensor<REAL> &endStress, int steps)
{
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> PrepareInitialMat *** Preconsolidating the material\n";
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
    REAL multipl;
    int i;
    TPZTensor<REAL> strain, localLoad, diffStress;
    
    diffStress = endStress;
    diffStress.Add(initialStress, -1.);
    
    for(i = 0; i <= steps; i++)
    {
        cout << "Starting step " << i << " of " << steps << endl;
        if(i == 0)
        {
            multipl = 0;
        }else
        {
            multipl = (REAL)i / (REAL)steps;
        }
        
        localLoad = initialStress;
        localLoad.Add(diffStress, multipl);
        
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "*** PrepareInitialMat *** Applying the " << i << "-th step with load = " << localLoad;
            sout << "\nPlasticState = " << mat.GetState();
            //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
        }
#endif
        
        mat.ApplyLoad(localLoad, strain);
        
#ifdef LOG4CXX
        {
            std::stringstream sout;
            sout << "*** PrepareInitialMat *** Achieved strain = " << strain;
            //		LOGPZ_INFO(PEPLogger,sout.str().c_str());
        }
#endif
        
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "<<< PrepareInitialMat *** End of material preconsolidation";
        sout << "\nMaterial Properties:\n";
        mat.Print(sout);
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
}

#include "pzstring.h"


REAL propPG(REAL ratio, int nrad, int ncirc, int i)
{	
    REAL q = 1.5;
    
    //REAL a0 = (ratio - 1.) * (q - 1.)/(pow(q, nrad+1) -1.);
    
    REAL SPi = (pow(q,i+1) - 1.) / (q - 1.) -1.;
    REAL SPn = (pow(q,nrad+1) - 1.) / (q - 1.) -1.;
    
    return SPi / SPn;
}


void QuarterWellboreGeom(int ncirc, 
                         REAL ioratio, 
                         TPZVec< TPZVec<REAL> > &pt, 
                         TPZVec< TPZVec<int> > &el, 
                         TPZVec< MElementType > &elType,
                         TPZVec< TPZString > &elName,
                         int & nrad)
{
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> QuarterWellboreGeom ***";
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
    REAL pi = 3.1415927;
    int i, j;
    ncirc = 2 * (ncirc / 2); //ensuring it's even
    int ncircpt = ncirc + 1;
    REAL q = 1.4; // 1.5
    nrad = static_cast<int>(( log(2.*((REAL) ncirc)*(ioratio - 1.)/pi)/log(q) -1. )*1.4 ); //1.6
    
    int nradpt = nrad + 1;
    int nel3d = ncirc * nrad;
    
    int nel = 3*nel3d + 2*ncirc + 2*nrad;
    
    ///nel = nel3d;////
    
    int nlayerpt = ncircpt * nradpt;
    
    pt.Resize(2 * nlayerpt);
    el.Resize(nel);
    elType.Resize(nel);
    elName.Resize(nel);
    
    //creating the vertices
    // wellbore vertices
    for( i = 0 ; i < ncircpt; i++)
    {
        REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
        TPZVec< REAL > pti(3,0.);
        pti[0] = cos(theta);
        pti[1] = sin(theta);
        pt[i] = pti;
    }
    // external x = ratio points
    for( i = 0; i <= ncirc / 2; i ++)
    {
        REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
        //REAL scale = ioratio / cos(theta);
        TPZVec< REAL > pti(3,0.);
        //pti[0] = ioratio;
        //pti[1] = sin(theta) * scale;
        pti[0] = ioratio * cos(theta);
        pti[1] = ioratio * sin(theta);
        pt[i + ncircpt * nrad] = pti;
    }
    // external y = ratio points
    for( i = ncirc / 2; i < ncircpt ; i ++)
    {
        REAL theta = pi/2.*((REAL) i ) / ((REAL) ncirc);
        //REAL scale = ioratio / sin(theta);
        TPZVec< REAL > pti(3,0.);
        //pti[0] = cos(theta) * scale;
        //pti[1] = ioratio;
        pti[0] = ioratio * cos(theta);
        pti[1] = ioratio * sin(theta);
        pt[i + ncircpt * nrad] = pti;
    }
    // filling the inbetween points
    for( i = 1; i < nrad; i++)
    {
        REAL prop = propPG(ioratio, nrad, ncirc, i);
        for( j = 0; j < ncircpt; j++)
        {
            TPZVec< REAL > pti(3,0.);
            TPZVec< REAL > ptExt(pt[ncircpt * nrad + j]);
            REAL sizeExt = ptExt[0] * ptExt[0] + ptExt[1] * ptExt[1] + ptExt[2] * ptExt[2];
            sizeExt = sqrt(sizeExt);
            pti[0] = pt[j][0] * (1. - prop) + ptExt[0] * prop * sqrt(ioratio / sizeExt);
            pti[1] = pt[j][1] * (1. - prop) + ptExt[1] * prop * sqrt(ioratio / sizeExt);
            pti[2] = pt[j][2] * (1. - prop) + ptExt[2] * prop * sqrt(ioratio / sizeExt);
            pt[i*ncircpt + j] = pti;
        }
    }
    // creating the z-layer points
    for( i = 0 ; i < nlayerpt; i++)
    {
        TPZVec< REAL > pti(3,0.);
        pti[0] = pt[i][0];
        pti[1] = pt[i][1];
        pti[2] = 1.;
        pt[i + nlayerpt] = pti;
    }
    //creating the 3D elements
    for(i = 0; i < nrad; i++)
    {
        for(j = 0; j < ncirc; j++)
        {
            int index = i * ncirc + j;
            TPZVec< int > eli(8,0);
            eli[0] = ncircpt * i + j;
            eli[1] = ncircpt * (i + 1) + j;
            eli[2] = eli[1] + 1;
            eli[3] = eli[0] + 1;
            eli[4] = eli[0] + nlayerpt;
            eli[5] = eli[1] + nlayerpt;
            eli[6] = eli[2] + nlayerpt;
            eli[7] = eli[3] + nlayerpt;
            el[index] = eli;
            elType[index] = ECube;
            elName[index] = "Interior";
        }
    }
    int lastIndex = nel3d;
    
    // Wellbore Faces
    for(i = 0; i < ncirc; i++)
    {
        TPZVec< int > eli(4);
        int index = i + lastIndex;
        
        eli[0] = el[i][0];
        eli[1] = el[i][3];
        eli[2] = el[i][7];
        eli[3] = el[i][4];
        
        el[index] = eli;
        elType[index] = EQuadrilateral;
        elName[index] = "Wellbore";
    }
    lastIndex += ncirc;
    
    //creating the 2D border elements
    for(i = 0; i < nel3d; i++)
    {
        TPZVec< int > elTop(4), elBot(4);
        int indexBot = i + lastIndex;
        int indexTop = i + lastIndex + nel3d;
        
        for(j = 0 ; j < 4; j++)
        {
            elBot[j] = el[i][j];
            elTop[j] = el[i][j+4];
        }
        el[indexBot] = elBot;
        elType[indexBot] = EQuadrilateral;
        elName[indexBot] = "Z- Plane Strain";
        
        el[indexTop] = elTop;
        elType[indexTop] = EQuadrilateral;
        elName[indexTop] = "Z+ Plane Strain";
    }
    lastIndex += 2* nel3d;
    
    // Lower Symmetry
    for(i = 0; i < nrad; i++)
    {
        TPZVec< int > eli(4);
        int index = i + lastIndex;
        int refIndex = ncirc*i;
        
        eli[0] = el[refIndex][0];
        eli[1] = el[refIndex][1];
        eli[2] = el[refIndex][5];
        eli[3] = el[refIndex][4];
        
        el[index] = eli;
        elType[index] = EQuadrilateral;
        elName[index] = "Lower Symmetry";
    }
    lastIndex += nrad;	
    
    // Left Symmetry
    for(i = 0; i < nrad; i++)
    {
        TPZVec< int > eli(4);
        int index = i + lastIndex;
        int refIndex = ncirc*(i+1) - 1;
        
        eli[0] = el[refIndex][3];
        eli[1] = el[refIndex][2];
        eli[2] = el[refIndex][6];
        eli[3] = el[refIndex][7];
        
        el[index] = eli;
        elType[index] = EQuadrilateral;
        elName[index] = "Left Symmetry";
    }
    lastIndex += nrad;
    
    // Right Farfield
    for(i = 0; i < ncirc / 2; i++)
    {
        TPZVec< int > eli(4);
        int index = i + lastIndex;
        int refIndex = (nrad-1)*ncirc + i;
        
        eli[0] = el[refIndex][1];
        eli[1] = el[refIndex][2];
        eli[2] = el[refIndex][6];
        eli[3] = el[refIndex][5];
        
        el[index] = eli;
        elType[index] = EQuadrilateral;
        elName[index] = "Right Farfield";
    }
    //lastIndex += ncirc / 2;
    
    // Top Farfield
    for(i = ncirc / 2; i < ncirc; i++)
    {
        TPZVec< int > eli(4);
        int index = i + lastIndex;
        int refIndex = (nrad-1)*ncirc + i;
        
        eli[0] = el[refIndex][1];
        eli[1] = el[refIndex][2];
        eli[2] = el[refIndex][6];
        eli[3] = el[refIndex][5];
        
        el[index] = eli;
        elType[index] = EQuadrilateral;
        elName[index] = "Top Farfield";
    }
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "<<< QuarterWellboreGeom ***";
        
        for(i = 0; i < 2 * nlayerpt; i++)
            sout << endl << "\tnode " << i << ": " << pt[i][0] << "\t" << pt[i][1] << "\t" << pt[i][2];
        for(i = 0; i < nel; i++)
        {
            sout << endl << "\tel " << i << ": " << elName[i].Str() << ";\tnodes:\t";
            int n = el[i].NElements();
            for(j = 0; j < n; j++)sout << el[i][j] << "\t";
        }
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
}


TPZCompMesh * CreateQuarterWellboreMesh( int gOrder,
                                        int ncirc,
                                        REAL ioratio,
                                        TPZMaterial * pMat,
                                        TPZFMatrix<REAL> & BCStressState,
                                        TPZFMatrix<REAL> & WellboreStressState,
                                        int allNeumannBC = 0)
{
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << ">>> CreateOneElCMesh ***";
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
    int matId = pMat->Id();
    int nstate = pMat->NStateVariables();
    int i;
    int nrad;
    //int ncirc = 4, nrad;
    //REAL ioratio = 10.;
    
    TPZVec< TPZVec< REAL > > pt;
    TPZVec< TPZVec< int > > el;
    TPZVec< MElementType > elType;
    TPZVec< TPZString > elName;
    TPZFMatrix<REAL> val1(nstate,nstate), val2(nstate,1);
    TPZFMatrix<REAL> bcNormal(nstate, 1, 0.);
    TPZFMatrix<REAL> bcStressState(nstate, nstate, 0.),
    wellboreStressState(nstate, nstate, 0.);
    
    QuarterWellboreGeom(ncirc, ioratio, pt, el, elType, elName, nrad);
    int nel = el.NElements();
    TPZGeoMesh * pGMesh = new TPZGeoMesh;
    
    // preparing nodes and elements
    int npt = pt.NElements();
    pGMesh->NodeVec().Resize(npt);
    for(i = 0; i < npt; i++)
        pGMesh->NodeVec()[i].Initialize(pt[i], *pGMesh);
    
    for(i = 0; i < nel; i++)
    {		
        matId = 0;
        if(!strcmp(elName[i].Str(),"Interior"))       matId =  1;
        if(!strcmp(elName[i].Str(),"Z- Plane Strain"))matId = -1;
        if(!strcmp(elName[i].Str(),"Z+ Plane Strain"))matId = -2;
        if(!strcmp(elName[i].Str(),"Lower Symmetry")) matId = -3;
        if(!strcmp(elName[i].Str(),"Left Symmetry"))  matId = -4;
        if(!strcmp(elName[i].Str(),"Right Farfield")) matId = -5;
        if(!strcmp(elName[i].Str(),"Top Farfield"))   matId = -5;
        if(!strcmp(elName[i].Str(),"Wellbore"))		  matId = -7;
        
        if(matId != 0)
        {
            pGMesh->CreateGeoElement(elType[i], el[i], matId, i);
        }else{
            PZError << "\nQuarterWellboreMesh error - element " << i << " without material assignement.\n";
        }
    }
    
    pGMesh->BuildConnectivity();
    
    TPZCompEl::SetgOrder(gOrder);
    
    // Creating Computation Mesh
    TPZCompMesh * pCMesh = new TPZCompMesh(pGMesh);
    
    TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(pCMesh);
    
    //Preparing Material
    pMat -> SetForcingFunction(NULL);
    TPZAutoPointer<TPZMaterial> mat(pMat);
    pCMesh->InsertMaterialObject(mat);
    
    //preparing BCs
    
    if(allNeumannBC)
    {
        for(i = -1; i > -7; i--)
        {
            TPZAutoPointer<TPZMaterial> bc;
            bc = mat->CreateBC(mat, i, 4 /*StressField Neumann BCType*/, BCStressState, val2);	
            pCMesh->InsertMaterialObject(bc);
        }
    }else
    {
        {//"Z- Plane Strain" bc -1
            bcNormal.Zero();
            bcNormal(2,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc;
            val2(2,0) = 1.;
            bc = mat->CreateBC(mat, -1, 3 /*Directional Dirichlet BCType*/, val1, val2);
            pCMesh->InsertMaterialObject(bc);
        }
        {//"Z+ Plane Strain" bc -2
            bcNormal.Zero();
            bcNormal(2,0) = 1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc;
            val2(2,0) = 1.;
            bc = mat->CreateBC(mat, -2, 3 /*Directional Dirichlet BCType*/, val1, val2);	
            pCMesh->InsertMaterialObject(bc);
        }
        {//"Lower Symmetry"  bc -3
            bcNormal.Zero();
            bcNormal(1,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc;
            val2(1,0) = 1.;
            bc = mat->CreateBC(mat, -3, 3 /*Directional Dirichlet BCType*/, val1, val2);	
            pCMesh->InsertMaterialObject(bc);
        }	
        {//"Left Symmetry"   bc -4
            bcNormal.Zero();
            bcNormal(0,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc;
            val2(0,0) = 1.;
            bc = mat->CreateBC(mat, -4, 3 /*Directional Dirichlet BCType*/, val1, val2);	
            pCMesh->InsertMaterialObject(bc);
        }
        {//"Right & Top Farfield"  bc -5
            val2.Zero();
            val2(0,0) = 1.;
            val2(1,0) = 1.;
            TPZAutoPointer<TPZMaterial> bc;
            bc = mat->CreateBC(mat, -5, 3, val1, val2);	//Directional Dirichlet BCType
            pCMesh->InsertMaterialObject(bc);
        }
        
    }
    
    {//"Wellbore" bc -7
        val2.Zero();
        TPZAutoPointer<TPZMaterial> bc;
        bc = mat->CreateBC(mat, -7, 4 , WellboreStressState, val2);//StressField Neumann
        pCMesh->InsertMaterialObject(bc);	
    }
    
    // building mesh connections
    
    
    //
    //preparing BCs
    if(allNeumannBC)
    {
        for(i = -1; i > -7; i--)
        {
            TPZAutoPointer<TPZMaterial> bc0;
            bc0 = mat->CreateBC(mat, i, 4 /*StressField Neumann BCType*/, BCStressState, val2);	;
            pCMesh->InsertMaterialObject(bc0);
        }
    }else
    {
        {//"Z- Plane Strain" bc -1
            bcNormal.Zero();
            bcNormal(2,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc1;
            val2(2,0) = 1.;
            bc1 = mat->CreateBC(mat, -1, 3 /*Directional Dirichlet BCType*/, val1, val2);
            pCMesh->InsertMaterialObject(bc1);
        }
        {//"Z+ Plane Strain" bc -2
            bcNormal.Zero();
            bcNormal(2,0) = 1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc2;
            val2(2,0) = 1.;
            bc2 = mat->CreateBC(mat, -2, 3 /*Directional Dirichlet BCType*/, val1, val2);
            pCMesh->InsertMaterialObject(bc2);
        }
        {//"Lower Symmetry"  bc -3
            bcNormal.Zero();
            bcNormal(1,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc3;
            val2(1,0) = 1.;
            bc3 = mat->CreateBC(mat, -3, 3 /*Directional Dirichlet BCType*/, val1, val2);	
            pCMesh->InsertMaterialObject(bc3);
        }	
        {//"Left Symmetry"   bc -4
            bcNormal.Zero();
            bcNormal(0,0) = -1.;
            val2.Zero();
            TPZAutoPointer<TPZMaterial> bc4;
            val2(0,0) = 1.;
            bc4 = mat->CreateBC(mat, -4, 3 /*Directional Dirichlet BCType*/, val1, val2);	
            pCMesh->InsertMaterialObject(bc4);
        }
        {//"Right & Top Farfield"  bc -5
            val2.Zero();
            
            val2(0,0) = 1.;
            val2(1,0) = 1.;
            TPZAutoPointer<TPZMaterial> bc5;
            bc5 = mat->CreateBC(mat, -5, 3, val1, val2);	//Directional Dirichlet BCType
            pCMesh->InsertMaterialObject(bc5);
        }
        
    }
    
    {//"Wellbore" bc -7
        val2.Zero();
        TPZAutoPointer<TPZMaterial> bc6;
        bc6 = mat->CreateBC(mat, -7, 4, WellboreStressState, val2);//StressField Neumann
        pCMesh->InsertMaterialObject(bc6);	
    }
    
    // building mesh connections
    
    pCMesh->AutoBuild();
    ofstream arg("MALHAPOCO.txt");
    pCMesh->Print(arg);
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        sout << "<<< CreateQuarterWellboreMesh *** with ncirc = " << ncirc << " * nrad = " << nrad << " elements";
        //	LOGPZ_INFO(PEPLogger,sout.str().c_str());
    }
#endif
    
    return pCMesh;

}


void BrazilianPlasticAnalysis2D()
{
	TPZFMatrix<REAL> BeginStress(3,3,0.), EndStress(3,3,0.), EndStress2(3,3,0.);
	TPZFMatrix<REAL> val1(3,1,0.);TPZFMatrix<REAL> val2(3,1,0.);TPZFMatrix<REAL> BeginForce(3,1,0.);TPZFMatrix<REAL> EndForce(3,1,0.);
	
	int BC1,BC2,nsteps,taxa,nnewton;
	int h=3;
	int order = 2;
	REAL tol = 1.e-5;
	nnewton = 3;
	//BC1 = -2;//down line
	//BC2 = -3;//upper line
    BC1 = -6;//down line
	BC2 = -5;//upper line
    
    

 // BC1 = -2;//EM CIMA -
//	BC2 = -1;//EM BAIXO +
	nsteps = 10;
	taxa = 1;
	BeginForce(1,0) =10.;
	EndForce(1,0) = 30.;
	
	
	TPZGeoMesh * MESH = new TPZGeoMesh;
	//void BrazilianTestGeoMesh::TransformBlendToLinearMesh2(TPZGeoMesh *newlinearmesh, int h)
	//BrazilianTestGeoMesh::TransformBlendToLinearMesh2(MESH,h);
	int refdir = 5;
	//MESH = BrazilianTestGeoMesh::TwoDMeshII(h,refdir);
    MESH = BrazilianTestGeoMesh::MalhaPredio();
	
	ofstream arg1("GeoMesh.txt");
	MESH->Print(arg1);
	TPZCompEl::SetgOrder(order);
	TPZCompMesh *CMESH = new TPZCompMesh(MESH);
	
    
	TPZElastoPlasticAnalysis::SetAllCreateFunctionsWithMem(CMESH);
    
    TPZLadeKim LK;
    TPZLadeKim::PlainConcrete(LK);
    
    TPZVonMises VM;
    TPZVonMises::Steel(VM);
    
    TPZDruckerPrager DP;
    TPZDruckerPrager::PlainConcreteMPa(DP);
    
    TPZDruckerPrager Very_RigidDP;
    TPZDruckerPrager::VeryRigidMaterial(Very_RigidDP);
    
    


	
    //TPZElasticityMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, int plainstress = 1);
    //TPZElasticityMaterial  matelast(2, 10000000., 0.1, 0.,0.,0);  

	TPZMatElastoPlastic2D<TPZDruckerPrager> EPMat1(1,1);
	EPMat1.SetPlasticity(DP);
    
    TPZMatElastoPlastic2D<TPZDruckerPrager> EPMat2(2,1);
	EPMat2.SetPlasticity(Very_RigidDP);
    
    
//    TPZMatElastoPlastic2D<TPZLadeKim> EPMat1(1,1);
//	EPMat1.SetPlasticity(LK);
    
 //   TPZMatElastoPlastic2D<TPZLadeKim> EPMat2(2,1);
//	EPMat2.SetPlasticity(LK);
    
    //TPZMatElastoPlastic2D<TPZVonMises> EPMat2(1,1);
	//EPMat2.SetPlasticity(VM);
    //TPZMatElastoPlastic2D<TPZLadeKim> EPMatx(2,1);
    //EPMatx.SetPlasticity(LK);
    
    
	TPZAutoPointer<TPZMaterial> plastic(&EPMat1);
    TPZAutoPointer<TPZMaterial> elastic(&EPMat2);
    
	
	//plastic->Print(cout);
	CMESH->InsertMaterialObject(plastic);
    CMESH->InsertMaterialObject(elastic);

    CMeshTwoMaterials(CMESH,elastic,plastic);

	
	ofstream arg("CMESHPLASTIC2D.txt");
	CMESH->Print(arg);
	
	TPZElastoPlasticAnalysis EPAnalysis(CMESH,cout);
	
	SolveSistII(EPAnalysis,CMESH);
    //EPAnalysis.SetBiCGStab(5000, 1.e-10);
   // EPAnalysis.SetBiCGStab_Jacobi(5000,1.e-12);
	
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
	TPZFStructMatrix structmatrix(PPAnalysis.Mesh());
	PPAnalysis.SetStructuralMatrix(structmatrix);
    //TPZVec<int> PostProcMatIds(1,1);
    TPZVec<int> PostProcMatIds(2);
    PostProcMatIds[0]=1;
    PostProcMatIds[1]=2;
	TPZVec<std::string> PostProcVars, scalNames, vecNames;
	SetUPPostProcessVariablesII(PostProcVars,scalNames,vecNames);
	PPAnalysis.SetPostProcessVariables(PostProcMatIds, PostProcVars);
	
	EPAnalysis.TransferSolution(PPAnalysis);
	
	cout << "\nDefining Graph Mesh\n";
	int dimension =2;
	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"BRLadep5p2.vtk");
	
	cout << "\nExporting First Solution without any refinement - initial solution might be smooth enough and a real mesh size output is of interest\n";
	
	
	PPAnalysis.PostProcess(0/*pOrder*/);
	
    //int pOrder =0;
	ManageIterativeProcessII(EPAnalysis,cout,tol,nnewton,BC1,BC2,nsteps,taxa,BeginStress,EndStress,BeginForce,EndForce,&PPAnalysis,0);

	//	cout << "\nInitial Solution Exported. Solving Problem\n";
	//	EPAnalysis.IterativeProcess(cout, 1.e-5, 30);
	//	cout << " \n PPPLASTIC-SOL \n";
	//	cout << EPAnalysis.Solution() << endl;
	//	EPAnalysis.AcceptSolution();
	//	EPAnalysis.TransferSolution(PPAnalysis);
	//	PPAnalysis.PostProcess(0);
	//	
	PPAnalysis.DefineGraphMesh(dimension,scalNames,vecNames,"BRLadep5p2.vtk");
	PPAnalysis.PostProcess(0);
	PPAnalysis.CloseGraphMesh();
	
}

void CMesh2DBCII(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
	
	
	int neumannPressure = 5;
	int neumann= 1;
	int dirichelet = 0;
	int mixed =2;
	REAL BIG = 1.e-12;
	
	
	//down line
	TPZFMatrix<REAL> k1(2,2,0.);
	TPZFMatrix<REAL> f1(2,1,0.);
	f1(1,0) = 17.;
	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -2, neumann, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
	
	//upper line
	TPZFMatrix<REAL> k2(2,2,0.);
	TPZFMatrix<REAL> f2(2,1,0.);
	f2(1,0) = -17.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -3, neumann, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
    
    
	TPZFMatrix<REAL> k6(2,2,0.);
	TPZFMatrix<REAL> f6(2,1,0.);
	TPZAutoPointer<TPZMaterial> ContBC6 = mat->CreateBC(mat, -8, 0, k6, f6);
	CMESH->InsertMaterialObject(ContBC6);
    
    TPZFMatrix<REAL> k7(2,2,0.);
	TPZFMatrix<REAL> f7(2,1,0.);
    //f7(1,0)=1.;
    k7(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC7 = mat->CreateBC(mat, -9, 2, k7, f7);
	CMESH->InsertMaterialObject(ContBC7);
    
    TPZFMatrix<REAL> k8(2,2,0.);
	TPZFMatrix<REAL> f8(2,1,0.);
    //f8(1,0)=1.;
    k8(1,1)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC8 = mat->CreateBC(mat, -10, 2, k8, f8);
	CMESH->InsertMaterialObject(ContBC8);
    
    TPZFMatrix<REAL> k3(2,2,0.);
	TPZFMatrix<REAL> f3(2,1,0.);
    //f3(0,0)=1.;
    k3(0,0)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -11, 2, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
	
    
	TPZFMatrix<REAL> k4(2,2,0.);
	TPZFMatrix<REAL> f4(2,1,0.);
    //f4(0,0)=1.;
    k4(0,0)=BIG;
	TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -12, 2, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
/*	

	TPZFMatrix<REAL> k3(2,2,0.);
	TPZFMatrix<REAL> f3(2,1,0.);
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -4, mixed, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
	
	TPZFMatrix<REAL> k4(3,3,0.);
	TPZFMatrix<REAL> f4(3,1,0.);
	TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -5, mixed, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
    
	TPZFMatrix<REAL> k5(3,3,0.);
	TPZFMatrix<REAL> f5(3,1,0.);
	TPZAutoPointer<TPZMaterial> ContBC5 = mat->CreateBC(mat, -6, mixed, k5, f5);
	CMESH->InsertMaterialObject(ContBC5);
    
	TPZFMatrix<REAL> k6(3,3,0.);
	TPZFMatrix<REAL> f6(3,1,0.);
	TPZAutoPointer<TPZMaterial> ContBC6 = mat->CreateBC(mat, -7, mixed, k6, f6);
	CMESH->InsertMaterialObject(ContBC6);
*/	
	CMESH->AutoBuild();
}

void CMeshTwoMaterials(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> matelast,TPZAutoPointer<TPZMaterial> matplastic)
{
	
    //Deslocamento na linha superior sinal negativo
    TPZFMatrix<REAL> k1(2,2,0.);
	TPZFMatrix<REAL> f1(2,1,0.);
    f1(1,0)=-1.;
	TPZAutoPointer<TPZMaterial> ContBC1 = matelast->CreateBC(matelast, -5, 1, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
    
    
    //Deslocamento na linha inferior sinal positvo
    TPZFMatrix<REAL> k2(2,2,0.);
	TPZFMatrix<REAL> f2(2,1,0.);
    f2(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC2 = matelast->CreateBC(matelast, -6, 1, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
    
    
    //Nos centrail
    TPZFMatrix<REAL> k3(2,2,0.);
	TPZFMatrix<REAL> f3(2,1,0.);
    f3(0,0) = 1.;
    f3(1,0) = 1.;
	TPZAutoPointer<TPZMaterial> ContBC3 = matplastic->CreateBC(matplastic, -1, 3, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
    
    //Nos direita
    TPZFMatrix<REAL> k4(2,2,0.);
	TPZFMatrix<REAL> f4(2,1,0.);
    f4(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC4 = matplastic->CreateBC(matplastic, -2, 3, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
    
    //Nos esquerda
    TPZFMatrix<REAL> k5(2,2,0.);
	TPZFMatrix<REAL> f5(2,1,0.);
    f5(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC5 = matplastic->CreateBC(matplastic, -3, 3, k5, f5);
	CMESH->InsertMaterialObject(ContBC5);
    
    
//    std::set<int> matids;
//	matids.insert(2);
//	matids.insert(1);
    CMESH->AutoBuild();
}



void CMeshGid(TPZCompMesh *CMESH, TPZAutoPointer<TPZMaterial> mat)
{
    
    //Deslocamento na linha superior sinal negativo
    TPZFMatrix<REAL> k1(2,2,0.);
	TPZFMatrix<REAL> f1(2,1,0.);
    f1(1,0)=-1.;
	TPZAutoPointer<TPZMaterial> ContBC1 = mat->CreateBC(mat, -5, 1, k1, f1);
	CMESH->InsertMaterialObject(ContBC1);
    
    
    //Deslocamento na linha inferior sinal positvo
    TPZFMatrix<REAL> k2(2,2,0.);
	TPZFMatrix<REAL> f2(2,1,0.);
    f2(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC2 = mat->CreateBC(mat, -6, 1, k2, f2);
	CMESH->InsertMaterialObject(ContBC2);
    
    
    //Nos centrail
    TPZFMatrix<REAL> k3(2,2,0.);
	TPZFMatrix<REAL> f3(2,1,0.);
   // f3(0,0) = 1.;
    //f3(1,0) = 1.;
	TPZAutoPointer<TPZMaterial> ContBC3 = mat->CreateBC(mat, -1, 0, k3, f3);
	CMESH->InsertMaterialObject(ContBC3);
    
    //Nos direita
    TPZFMatrix<REAL> k4(2,2,0.);
	TPZFMatrix<REAL> f4(2,1,0.);
    //f4(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC4 = mat->CreateBC(mat, -2, 0, k4, f4);
	CMESH->InsertMaterialObject(ContBC4);
    
    //Nos esquerda
    TPZFMatrix<REAL> k5(2,2,0.);
	TPZFMatrix<REAL> f5(2,1,0.);
   // f5(1,0)=1.;
	TPZAutoPointer<TPZMaterial> ContBC5 = mat->CreateBC(mat, -3, 0, k5, f5);
	CMESH->InsertMaterialObject(ContBC5);
    
   
	CMESH->AutoBuild();
 
}


void SetUPPostProcessVariablesII(TPZVec<std::string> &postprocvars, TPZVec<std::string> &scalnames, TPZVec<std::string> &vecnames )
{
	
	scalnames.Resize(10);
	scalnames[0] = "Alpha";
	scalnames[1] = "PlasticSteps";
	scalnames[2] = "VolElasticStrain";
	scalnames[3] = "VolPlasticStrain";
	scalnames[4] = "VolTotalStrain";
	scalnames[5] = "I1Stress";
	scalnames[6] = "J2Stress";
	scalnames[7] = "YieldSurface";
	scalnames[8] = "EMisesStress";
	scalnames[9] = "TotalPlasticStrain";
    
  /*  
    scalnames.Resize(9);
	scalnames[0] = "Alpha";
	scalnames[1] = "PlasticSteps";
	scalnames[2] = "VolElasticStrain";
	scalnames[3] = "VolPlasticStrain";
	scalnames[4] = "VolTotalStrain";
	scalnames[5] = "I1Stress";
	scalnames[6] = "J2Stress";
	scalnames[7] = "EMisesStress";
	scalnames[8] = "TotalPlasticStrain";
	*/
    
	vecnames.Resize(5);
	vecnames[0] = "Displacement";
	vecnames[1] = "NormalStress";
	vecnames[2] = "ShearStress";
	vecnames[3] = "NormalStrain";
	vecnames[4] = "ShearStrain";
	//vecnames[5] = "NormalPlasticStrain";
    
	
    
//    if(!strcmp("displacement",name.c_str()))     return 9;
//	if(!strcmp("Pressure",name.c_str()))         return 1;
//	if(!strcmp("MaxStress",name.c_str()))        return 2;
//	if(!strcmp("PrincipalStress1",name.c_str())) return 3;
//	if(!strcmp("PrincipalStress2",name.c_str())) return 4;
//	if(!strcmp("SigmaX",name.c_str()))           return 5;
//	if(!strcmp("SigmaY",name.c_str()))           return 6;
//	if(!strcmp("TauXY",name.c_str()))            return 8;//Cedric
//	if(!strcmp("sig_x",name.c_str()))            return 5;
//	if(!strcmp("sig_y",name.c_str()))            return 6;
//	if(!strcmp("tau_xy",name.c_str()))           return 8;//Cedric
//	if(!strcmp("Displacement6",name.c_str()))    return 7;
//	if(!strcmp("Stress",name.c_str()))           return 10;
    
	
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


void ManageIterativeProcessII(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
							int BCId,int BCId2, int nsteps, REAL PGRatio,
							TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
							TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
							TPZPostProcAnalysis * ppAnalysis, int res)
{
	if(!analysis.Mesh())return;
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "<<< TPZElastoPlasticAnalysis::ManageIterativeProcess() ***";
		sout << "\nWith parameters:\n";
		sout << "\ntol = " << tol;
		sout << "\nnumiter = " << numiter;
		sout << "\nBCId = " << BCId;
		sout << "\nBCId2 = " << BCId2;
		sout << "\nnsteps = " << nsteps;
		sout << "\nPGRatio = " << PGRatio;
		sout << "\nval1Begin = " << val1Begin;
		sout << "\nval1End = " << val1End;
		sout << "\nval2Begin = " << val2Begin;
		sout << "\nval2End = " << val2End;
		if(ppAnalysis)
		{
			sout << "\nppanalysis set";
		}else
		{
			sout << "\nppanalysis NOT set";
		}
		//LOGPZ_INFO(logger,sout.str().c_str());
	}
#endif
	
	// computing the initial value for the PG progression such that its sum equals one;
	REAL a0;
	
	if(fabs(PGRatio - 1.) < 1.e-3)
	{
	    a0 = 1. / REAL(nsteps);
	}
	
	else
	{
		a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
	}
	TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
	TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
	
	deltaVal1 = val1End;
	deltaVal1.ZAXPY(-1., val1Begin);
	deltaVal2 = val2End;
	deltaVal2.ZAXPY(-1., val2Begin);	
	
	
	//-19
	TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
	if(!pBC)return;
	
	//-20
	//COMENTAR NO CILINDRO
		TPZAutoPointer<TPZMaterial> mat2 = analysis.Mesh()->FindMaterial(BCId2);
		TPZBndCond * pBC2 = dynamic_cast<TPZBndCond *>(mat2.operator->());
		if(!pBC2)return;
	
	bool linesearch= true;
    int i;
	for(i = 0; i < nsteps; i++)
	{
		REAL stepLen;
		if(fabs(PGRatio - 1.) < 1.e-3)
		{
			stepLen = REAL(i+1) / REAL(nsteps);
		}
		else
		{
		    stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
		}
		
		val1 = val1Begin;
		val1.ZAXPY(stepLen, deltaVal1);
		val2 = val2Begin;
		val2.ZAXPY(stepLen, deltaVal2);
		
		pBC->Val1() =1.* val1;
		pBC->Val2() =1.* val2;
		
     //   cout <<  "PRESSURE = " << val1 << endl;
		cout <<  "\nPRESSURE = " << val2 << endl;
		
		//COMENTAR NO CILINDRO
        pBC2->Val1() = -1.*val1;
        pBC2->Val2() = -1.*val2;		
		
		analysis.IterativeProcess(out, tol, numiter);
//        virtual void IterativeProcess(std::ostream &out,REAL tol,int numiter, bool linesearch = false, bool checkconv = false);
 

        if(i==0/*(nsteps/2)*/){
            linesearch = false;
        }
        
  //      analysis.IterativeProcess(out, tol,numiter,linesearch,false);
        
		analysis.AcceptSolution();
		
		
		if(ppAnalysis)
		{
			analysis.TransferSolution(*ppAnalysis);
			ppAnalysis->PostProcess(res);
		}
	}
}


void ManageIterativeProcessIII(TPZElastoPlasticAnalysis &analysis ,std::ostream &out,REAL tol,int numiter,
                                                      int BCId, int nsteps, REAL PGRatio,
                                                      TPZFMatrix<REAL> & val1Begin, TPZFMatrix<REAL> & val1End,
                                                      TPZFMatrix<REAL> & val2Begin, TPZFMatrix<REAL> & val2End,
                                                      TPZPostProcAnalysis * ppAnalysis, int res)
{
	// computing the initial value for the PG progression such that its sum equals one;
	REAL a0;
	if(fabs(PGRatio - 1.) < 1.e-3)
	{
	    a0 = 1. / REAL(nsteps);
	}else{
		a0 = (PGRatio - 1) / (pow(PGRatio,nsteps) - 1.);
	}
	TPZFNMatrix<36> val1(6,6,0.), deltaVal1(6,6,0.);
	TPZFNMatrix< 6> val2(6,1,0.), deltaVal2(6,1,0.);
	
	deltaVal1 = val1End;
	deltaVal1.ZAXPY(-1., val1Begin);
	deltaVal2 = val2End;
	deltaVal2.ZAXPY(-1., val2Begin);
	
	// ZAXPY operation: *this += alpha * p			
    
	TPZAutoPointer<TPZMaterial> mat = analysis.Mesh()->FindMaterial(BCId);
	TPZBndCond * pBC = dynamic_cast<TPZBndCond *>(mat.operator->());
	if(!pBC)return;
    
    int i;
	for(i = 0; i < nsteps; i++)
	{
		REAL stepLen;
		if(fabs(PGRatio - 1.) < 1.e-3)
		{
			stepLen = REAL(i+1) / REAL(nsteps);
		}else{
		    stepLen = a0 * (pow(PGRatio,i+1) - 1) / (PGRatio - 1.);
		}
		
		val1 = val1Begin;
		val1.ZAXPY(stepLen, deltaVal1);
		val2 = val2Begin;
		val2.ZAXPY(stepLen, deltaVal2);
		
		pBC->Val1() = val1;
		pBC->Val2() = val2;
    
		
		analysis.IterativeProcess(out, tol, numiter);

		
		analysis.AcceptSolution();
		
		if(ppAnalysis)
		{

			analysis.TransferSolution(*ppAnalysis);
			ppAnalysis->PostProcess(res);
		}
	}
    

}



void SolveSistII(TPZAnalysis &an, TPZCompMesh *fCmesh)
{

    //TPZFStructMatrix full(fCmesh)
	TPZSkylineStructMatrix full(fCmesh);
	an.SetStructuralMatrix(full);
    
    
	TPZStepSolver<REAL> step;
  //  step.SetDirect(ELDLt);
  //  step.SetJacobi(5000, 1.e-12,0);
    step.SetDirect(ECholesky);
    
	an.SetSolver(step);
	
}