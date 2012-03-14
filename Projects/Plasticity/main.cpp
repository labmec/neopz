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


void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
                            int BCId,int BCId2, int nsteps, REAL PGRatio,
                            TPZFMatrix & val1Begin, TPZFMatrix & val1End,
                            TPZFMatrix & val2Begin, TPZFMatrix & val2End)
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
    TPZStepSolver step;
    full.SetNumThreads(8);
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
    
    
    TPZCompMesh * pCMesh = CreateQuarterWellboreMesh(pOrder, ncirc, ioRatio, &EPMat, BeginStress, EndStress, 0);
    
    //building analysis
    TPZElastoPlasticAnalysis EPAnalysis(pCMesh, std::cout);
    
    SolveSistLin2(EPAnalysis,pCMesh);
    EPAnalysis.SetBiCGStab(5000, 1.e-12);
    
    //    void ManageIterativeProcess(TPZElastoPlasticAnalysis &analysis , std::ostream &out,REAL tol,int numiter,
    //                                int BCId,int BCId2, int nsteps, REAL PGRatio,
    //                                TPZFMatrix & val1Begin, TPZFMatrix & val1End,
    //                                TPZFMatrix & val2Begin, TPZFMatrix & val2End)
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
        TPZFMatrix val2(4,1,0.), val1(4,4,0.);
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
    int choice;
    cin >> choice;
    
    switch(choice)
    {
        case(0):
            cout << "\n Choose the Plastic model tou need to run Isotropic compression: ";
            cout << "\n0 - Lade - Kim ";
            cout << "\n1 - Sandler Dimaggio ";
            cout << "\n2 - Drucker Prager ";
            cout << "\n";
            int choice2;
            cin >> choice2;
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




int main()
{
    //InitializeLOG();
  //  LKIsotropicCompression();
  //  SandlerDimaggioIsotropicCompression();
  //   DruckerIsotropicCompression();
   // LKBiaxialTest();
    
    cout << "\nPlease enter test type:";
    cout << "\n0 Material Point test ";
    cout << "\n1 Finite Elements test ";
    
    int testnum;
   // testnum =3;
    cin >> testnum;
    
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
                    pDP->fYC.SetUp(/*phi*/ 29.7/180. *M_PI ,/*innerMCFit*/1);
                    pDP->fTFA.SetUp(/*yield- coesao inicial*/ 12.8, /*k Modulo de hardening da coesao equivante 10^-3 Mpa a cada 0.1% de deformacao */1.);
                    pDP->fER.SetUp(/*young*/ 29269., /*poisson*/ 0.203);
                    pMat = pDP;
                    fileName << "_PRDPInscMPa";
                    loadMultipl = -1/145.03773801;
                    break;
                case(10):
                    pDP = new TPZDruckerPrager();
                    pDP->fYC.SetUp(/*phi*/ 29.7/180. * M_PI ,/*innerMCFit*/0);
                    pDP->fTFA.SetUp(/*yield- coesao inicial*/ 12.8, /*k Modulo de hardening da coesao equivante 10^-3 Mpa a cada 0.1% de deformacao */1.);
                    pDP->fER.SetUp(/*young*/ 29269., /*poisson*/ 0.203);
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
    
    
    
    return EXIT_SUCCESS;
    
    
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
                                        TPZFMatrix & BCStressState,
                                        TPZFMatrix & WellboreStressState,
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
    TPZFMatrix val1(nstate,nstate), val2(nstate,1);
    TPZFMatrix bcNormal(nstate, 1, 0.);
    TPZFMatrix bcStressState(nstate, nstate, 0.),
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
