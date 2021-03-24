/**
 * @file
 * @brief Contains the implementation of the TPZTimeTemp methods. 
 */

#include "TPZTimeTemp.h"
#include "pzlog.h"
#include "pzstack.h"
#include "pzvec.h"

TPZTimeTemp tempo;

TPZTimeTemp::TPZTimeTemp()				
{
	fNumEq = 0;										
	fNumthreads = 0;
	fPolyOrder = 0;
	fNumSub = 0;
	fNumEqCoarse = 0;
	fniter = 0;
	fnMultiply = 0;
	fnPreCond = 0;
	ft0sub = 0;								
	ft1comput = 0;									
	ft2congraph = 0;								
	ft3analysegraph = 0;							
	ft4identcorner = 0;								
	ft5dohrassembly = 0;
	ft55decompmatriznosinternos = 0;
	ft6iter = 0;	
	
}


void TPZTimeTemp::PrintHeader(std::ostream &out)											
{
	out << "Time_to_SubMesh\tTime_to_ComputeSystemofEq\tTime_to_ConvertGraph\tTime_to_AnalyseGraph\tTotal_time_to_IdentConerNodes\tTime_threadDohrmanAssembly" <<
	"\tTime_to_Decompose_InnerNodesMatrix\tTime_for_Iterations\tNumEq\tNum_Eq_Coarse\tNumber_of_Elements\tNum_threads\tPoly_Order\tNumSub\tNum_of_Iterations\tElemnts_fnMultiply\tElements_fnPreCond\n";
}

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
	
	out << ft0sub <<  "\t" << ft1comput <<  "\t" << ft2congraph <<  "\t" << ft3analysegraph <<  "\t" << ft4identcorner <<  "\t"  << ft5dohrassembly  <<  "\t"  << ft55decompmatriznosinternos <<  "\t"  << ft6iter
	<<  "\t" << fNumEq << "\t" << fNumEqCoarse << "\t" << fNumberofElements << "\t"  << fNumthreads << "\t"  << fPolyOrder << "\t" << fNumSub << "\t"  << fniter << "\t" << fnMultiply << "\t" 
	<< fnPreCond << "\t" << fMultiply << "\t" <<  fPreCond << std::endl;
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
	
	ReadFile >> ft0sub >> ft1comput >> ft2congraph >> ft3analysegraph >> ft4identcorner >> ft5dohrassembly >> ft55decompmatriznosinternos >> ft6iter >> fNumEq >> fNumEqCoarse 
	>> fNumberofElements >> fNumthreads >> fPolyOrder >> fNumSub >> fniter >> fnMultiply >> fnPreCond;
	
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