
#ifdef HAVE_CONFIG_H
#include <config.h>
#endif


#include "pzskylstrmatrix.h"


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
#include "pzmat2dlin.h"
#include "pzporoanalysis.h"


#include "TPZTensor.h"

#include "BrazilianTestGeoMesh.h"

#include "pzelast3d.h"

#include "pzcompelpostproc.h"
#include "pzpostprocmat.h"
#include "pzpostprocanalysis.h"

#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "pzbdstrmatrix.h"
#include "pzstepsolver.h"
#include <sstream>




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

#include "pzelastoplastic2D.h"
#include "tpzycvonmisescombtresca.h"
#include "TPZMohrCoulombNeto.h"
#include "TPZSandlerDimaggio.h"


void VisualizeSandlerDimaggio(std::stringstream &FileName, TPZSandlerDimaggio *pSD);

using namespace pzshape; // needed for TPZShapeCube and related classes



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
 
    
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
    

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
    
	TPZPostProcAnalysis PPAnalysis(&EPAnalysis);
    
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

int main()
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
	cout << "\n2) Graphical representation of Yield surface";
    
	//cin >> testNumber;
    testNumber = 2;
	
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
    
//    cin >> matNumber;
    matNumber = 7;
	
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
			pDP->fYC.SetUp(/*phi*/ 29.7/180. * pi ,/*innerMCFit*/1);
			pDP->fTFA.SetUp(/*yield- coesao inicial*/ 12.8, /*k Modulo de hardening da coesao equivante 10^-3 Mpa a cada 0.1% de deformacao */1.);
			pDP->fER.SetUp(/*young*/ 29269., /*poisson*/ 0.203);
			pMat = pDP;
		    fileName << "_PRDPInscMPa";
		    loadMultipl = -1/145.03773801;
            break;
		case(10):
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
    plasticTol = 0.1;
    
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
        case 2:
            copyStr << fileName.str();
            fileName << "Yield" << copyStr.str();
            if (pSD) {
                VisualizeSandlerDimaggio(fileName,pSD);
            }
            break;
		default:
			cout << "\nUnhandled Test Type. Exiting...";
            delete pMat;
	}
	
    return EXIT_SUCCESS;
    
    
}

#include "TPZGenSpecialGrid.h"

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio *pSD);

void VisualizeSandlerDimaggio(std::stringstream &fileName, TPZSandlerDimaggio *pSD)
{
    TPZGeoMesh *gmesh = TPZGenSpecialGrid::GeneratePolygonalSphereFromOctahedron(1., 0.001);
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

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio *pSD, TPZVec<REAL> &stressresult);

void BuildPlasticSurface(TPZCompMesh *cmesh, TPZSandlerDimaggio *pSD)
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

int ComputeMultiplier(TPZVec<REAL> &stress, TPZSandlerDimaggio *pSD,TPZVec<REAL> &stressresult)
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
    SANDLERDIMAGGIOPARENT *pSDP = dynamic_cast<SANDLERDIMAGGIOPARENT *>(pSD);
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
