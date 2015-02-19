/*
 *  pznondarcyanalysis.cpp
 *  PZ
 *
 *  Created by Omar Duran Triana on 5/21/14.
 *  Copyright 2014 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZMultilayerDarcyAnalysis.h"


TPZDarcyAnalysis::TPZDarcyAnalysis(TPZCompMesh *Computationalmesh ,std::ostream &out) : TPZNonLinearAnalysis(Computationalmesh,out)
{
	if(!Computationalmesh) DebugStop();
    fResidualLastState.Resize(1, 1);
    fResidualLastState.Zero();
    fSimulationData = NULL;

}


TPZDarcyAnalysis::~TPZDarcyAnalysis()
{
	
}

void TPZDarcyAnalysis::SetLastState()
{

}

void TPZDarcyAnalysis::SetNextState()
{

}

void TPZDarcyAnalysis::Assemble()
{

}

void TPZDarcyAnalysis::AssembleLastStep()
{

}

void TPZDarcyAnalysis::AssembleResidual()
{

}

void TPZDarcyAnalysis::Residual(TPZFMatrix<STATE> &residual, int icase)
{
	TPZNonLinearAnalysis::Residual(residual, icase);
	residual = fResidualLastState + residual;
}

void TPZDarcyAnalysis::ComputeTangent(TPZFMatrix<STATE> &tangent, TPZVec<REAL> &coefs, int icase)
{
	this->SetNextState();
	TPZDarcyAnalysis::ComputeTangent(tangent, coefs, icase);
}

void TPZDarcyAnalysis::TimeForward(TPZFMatrix<STATE> &AlphasAtNplusOne, TPZFMatrix<STATE> &AlphasAtN)
{
	this->LoadSolution(AlphasAtN);
	this->AssembleLastStep();
	
	REAL tol = 1.e-6;
	int numiter = 2;
	
	this->IterativeProcess(std::cout,tol,numiter,false,false);
	
	AlphasAtNplusOne = fSolution;
}

void TPZDarcyAnalysis::InitializeFirstSolution(TPZFMatrix<STATE> &AlphasAtN, REAL &ReferencePressure)
{

}

void TPZDarcyAnalysis::Run()
{
	
}

