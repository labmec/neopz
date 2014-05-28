/*
 *  pznondarcyanalysis.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "pznondarcyanalysis.h"
#include "tpzautopointer.h"
#include "pzcmesh.h"
#include "pzintel.h"
#include "pzstepsolver.h"


TPZNonDarcyAnalysis::TPZNonDarcyAnalysis(TPZCompMesh *cmesh,                                             
                                             std::ostream &out)
: TPZNonLinearAnalysis(cmesh,out),
fResLastState(0.)
{
	if(!cmesh) DebugStop();
}


TPZNonDarcyAnalysis::~TPZNonDarcyAnalysis()
{
	
}

void TPZNonDarcyAnalysis::SetMaterialLastState()
{
	globData.fState = ELastState;
}

void TPZNonDarcyAnalysis::SetMaterialNextState()
{
	globData.fState	= ECurrentState;
}

void TPZNonDarcyAnalysis::Assemble()
{
	this->SetMaterialNextState();
	TPZAnalysis::Assemble();
	fRhs += fResLastState;
}

void TPZNonDarcyAnalysis::AssembleLastStep()
{
	this->SetMaterialLastState();
	TPZAnalysis::AssembleResidual();
	fResLastState = fRhs;
}

void TPZNonDarcyAnalysis::AssembleResidual()
{
	this->SetMaterialNextState();
	TPZAnalysis::AssembleResidual();
	fRhs += fResLastState;
}

void TPZNonDarcyAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase)
{
	TPZNonLinearAnalysis::Residual(residual, icase);
	residual += fResLastState;
}

void TPZNonDarcyAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase)
{
	this->SetMaterialNextState();
	TPZNonLinearAnalysis::ComputeTangent(tangent, coefs, icase);
}

void TPZNonDarcyAnalysis::TimeStep(TPZFMatrix<STATE> &PressureN, TPZFMatrix<STATE> &PressureNp1)
{
	this->LoadSolution(PressureN);
	this->AssembleLastStep();
	
	REAL tol = 1.e-6;
	int numiter = 2;
	
	this->IterativeProcess(std::cout,tol,numiter,false,false);
	
	PressureNp1 = fSolution;
}

void TPZNonDarcyAnalysis::InitializeFirstSolution(TPZFMatrix<STATE> &Pressure, REAL &pini)
{

	Pressure.Redim(fSolution.Rows(),fSolution.Cols());
	TPZBlock<STATE> &block = this->Mesh()->Block();
	int nel = this->Mesh()->NElements();
	for (int iel = 0 ; iel < nel ; iel++){
		TPZCompEl *cel = this->Mesh()->ElementVec()[iel];
		if (!cel) continue;
		TPZInterpolatedElement *intel = dynamic_cast <TPZInterpolatedElement *> (cel);
		if (!intel) DebugStop();
		if (intel->Reference()->Dimension() != 2) continue; //soh elem 2d.
		int ncc = intel->NCornerConnects();
		//if (ncc != 4) DebugStop(); // I expect only quad element, althought it would work for triang
		for (int icc = 0 ; icc < ncc ; icc++){
			TPZConnect *con = &intel->Connect(icc);
			int conseq = con->SequenceNumber();
			int pos = block.Position(conseq);
			Pressure(pos,0) = pini;
		}
	}
}

void TPZNonDarcyAnalysis::RunAll()
{
	
	TPZFMatrix<STATE> PressureN, PressureNp1;
	REAL t =  0;
	this->InitializeFirstSolution(PressureN,globData.fPini);
	PressureNp1 = PressureN;
	
	this->LoadSolution(PressureN);
	
	// PostProcVars
	TPZStack<std::string> scalnames, vecnames;
	scalnames.Push("Pressure");
	vecnames.Push("Velocity");
	const int dim = 2;
	const int postprocessresolution = 2;
	std::string file = "NonDarcyFlow.vtk";
	this->DefineGraphMesh(dim, scalnames, vecnames, file);
	this->PostProcess(postprocessresolution);

	const int nsteps = globData.fNSteps;
	
	for (int istep = 0 ; istep < nsteps ; istep++){
		const REAL dt = globData.fDeltaT;
		t += dt;
		
		this->TimeStep(PressureN,PressureNp1);
		this->PostProcess(postprocessresolution);
		
		PressureN = PressureNp1;
	}
}

