#ifndef TPZTIMETEMP_H
#define TPZTIMETEMP_H


/*
 *  TPZTimeTemp.h
 *  SubStruct
 *
 *  Created by Bandit on 7/27/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <string>
#include <iostream>
#include <fstream>
#include "pzstack.h"



/// Tomada de tempos
class TPZTimeTemp
{
public:

	/// Number of equations
	int fNumEq;
	/// Number of Threads to decompose
	int fNumthreads_decompose;	
	/// Number of Submeshs assembled at the same time
	int fnumthreads_nsubmesh_assemble;
	/// Number of Threads to assemble
	int fNumthreads_assemble;
	/// Number of threads to multiply
	int fNumthreads_multiply;
	/// polynomial order
	int fPolyOrder;
	/// Number of Substructures
	int fNumSub;											
	/// Numero de equacoes coarse
	int fNumEqCoarse;									
	/// Number of iterations
	int fniter;
	/// Number of elements in fMultiply
	int fnMultiply;
	/// Number of elements in fPreCond
	int fnPreCond;	
	/// Number of elements in the mesh
	int fNumberofElements;
	/// Time for Substructuring Mesh
	REAL ft0sub;										
	/// Time for Computing the system of equations for each substructure
	REAL ft1comput;									
	/// Time for Convert Graph
	//REAL ft2congraph;								
	/// Time for AnalyseGraph
	//REAL ft3analysegraph;							
	/// Total Time for Identifying Corner Nodes = Convert Graph + AnalyseGraph
	REAL ft4identcorner;								
	/// Time for ThreadDohrmanAssembly 
	REAL ft5dohrassembly;
	/// Time for ThreadDohrmanDecompose
	REAL ft51dohrdecompose;
	/// Time to decompose the inner nodes matrix
	//REAL ft55decompmatriznosinternos; // nao usado!
	/// Total Time for Iterations
	REAL ft6iter;
	
	/// Time to Multiply for each iteration
	TPZStack<REAL> fMultiply;		
	/// Time to PreCond for each iteration
	TPZStack<REAL> fPreCond;	
	
	/// inicializa os atributos
	TPZTimeTemp();												
	
	/*
	/// Imprime o arquivo de tempos com titulos
	void Print(std::ostream &out);	
	 */
	
	/// Append to the file a line with the informations 
	void PrintLine(std::ostream &out);					
	
	/// Read one line of the file
	void ReadLine(std::istream &in);
		
	/// Print the header
	static void PrintHeader(std::ostream &out);	
	
	/// If the file does not exists, returns "true"
	static bool NeedsHeader(std::string &FileName);
	
private:
	
};


extern TPZTimeTemp tempo;



#endif

