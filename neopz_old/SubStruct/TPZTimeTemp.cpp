/*
 *  TPZTimeTemp.cpp
 *  SubStruct
 *
 *  Created by Bandit on 7/27/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZTimeTemp.h"
#include "pzlog.h"
#include "pzstack.h"
#include "pzvec.h"

TPZTimeTemp tempo;

TPZTimeTemp::TPZTimeTemp()				
{
	fNumEq = 0;										
	fNumthreads_assemble = 0;
	fNumthreads_multiply = 0;
	fNumthreads_decompose = 0;
	fnumthreads_nsubmesh_assemble = 0;
	fPolyOrder = 0;
	fNumSub = 0;
	fNumEqCoarse = 0;
	fniter = 0;
	fnMultiply = 0;
	fnPreCond = 0;
	ft0sub = 0;								
	ft1comput = 0;									
	//ft2congraph = 0;								
	//ft3analysegraph = 0;							
	ft4identcorner = 0;								
	ft5dohrassembly = 0;
	ft51dohrdecompose = 0;
	//ft55decompmatriznosinternos = 0;
	ft6iter = 0;	

}


void TPZTimeTemp::PrintHeader(std::ostream &out)											
{
	out << "Time_to_SubMesh\tTime_to_ComputeSystemofEq\tTime_to_IdentConerNodes\tTime_threadDohrmanAssembly" <<
	"\tTime_threadDohrmnaDecompose\tTime_for_Iterations\tPoly_Order\tNumEq\tNum_Eq_Coarse\tNumber_of_Elements\tNumSub\tNum_Submeshs_assembled_same_time\tNum_threads_assemble_each_submesh" <<
	"\tNum_threads_decompose\tNum_threads_multiply\tNum_of_Iterations\tElemnts_fnMultiply\tElements_fnPreCond\n";
}


/*
void TPZTimeTemp::Print(std::ostream &out)		
{
	
	if (!out)
	{
		std::cout << __PRETTY_FUNCTION__ << "Arquivo nao valido" << std::endl;
		return;
	}
	out.precision(5);
	out << "Time for Substructuring Mesh: " << ft0sub
	<< "\nTime for Computing the system of equations for each substructure: " << ft1comput
	<< "\nInside Identifying Corner Nodes:\n" 
	<< "Time for Convert Graph: " << ft2congraph
	<< "\nTime for AnalyseGraph: " << ft3analysegraph
	<< "\nTotal Time for Identifying Corner Nodes: " << ft4identcorner
	<< "\nTime for ThreadDohrmanAssembly: " << ft5dohrassembly
	<< "\nTotal Time for Iterations: " << ft6iter
	<< "\nTime to Multiply for each iteration: " << fMultiply
	<< "\nTime to PreCond for each iteration: " << fPreCond
	<< "\nNumber of Equations: " << fNumEq
	<< "\nNumber of Coarse Equations: " << fNumEqCoarse
	<< "\nNumber of Threads: " << fNumthreads
	<< "\nNumber of Substructures: " << fNumSub << "\n" << std::endl;
}
*/


void TPZTimeTemp::PrintLine(std::ostream &out)									
{
	
	out.precision(5);
	 
	if (!out)
	{
		std::cout << __PRETTY_FUNCTION__ << "Arquivo nao valido" << std::endl;
		return;
	}
	fnMultiply = fMultiply.NElements();
	fnPreCond = fPreCond.NElements();
	fniter = fnPreCond;
	
	out << ft0sub <<  "\t" << ft1comput <<  "\t" << ft4identcorner <<  "\t"  << ft5dohrassembly  <<  "\t"  << ft51dohrdecompose <<  "\t"  << ft6iter 
	<< "\t" << fPolyOrder << "\t" << fNumEq << "\t" << fNumEqCoarse << "\t" << fNumberofElements << "\t" << fNumSub << "\t" << fnumthreads_nsubmesh_assemble 
	<< "\t" << fNumthreads_assemble << "\t" << fNumthreads_decompose << "\t" << fNumthreads_multiply << "\t" << fniter << "\t" << fnMultiply 
	<< "\t" << fnPreCond << "\t" << fMultiply << "\t" <<  fPreCond << std::endl;
}

bool TPZTimeTemp::NeedsHeader(std::string &FileName)
{
	std::ifstream Toto (FileName.c_str());
	bool shouldprintheader = false;
	if (!Toto)
	{
		char buf[1024];
		Toto >> buf;
		if(!Toto) shouldprintheader = true;
	}
	return shouldprintheader;
}

void TPZTimeTemp::ReadLine(std::istream &ReadFile)
{

	ReadFile >> ft0sub >> ft1comput >> ft4identcorner >> ft5dohrassembly >> ft51dohrdecompose >> ft6iter >> fPolyOrder >> fNumEq >> fNumEqCoarse 
	>> fNumberofElements >> fNumSub >> fnumthreads_nsubmesh_assemble >> fNumthreads_assemble >> fNumthreads_decompose >> fNumthreads_multiply >> fniter >> fnMultiply >> fnPreCond;
	
	REAL temp;
	for (int j = 0; j < fnMultiply; j++)				
	{
		ReadFile >> temp;
		fMultiply.Push(temp);
	}
	
	for (int k = 0; k < fnPreCond; k++) 
	{
		ReadFile >> temp;
		fPreCond.Push(temp);
	}	
}