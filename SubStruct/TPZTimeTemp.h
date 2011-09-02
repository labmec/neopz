/**
 * @file
 * @brief Contains the TPZTimeTemp class which takes times.
 */
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

/**
 * @brief Takes times. (Tomada de tempos). \ref util "Utility"
 * @ingroup util
 */
class TPZTimeTemp
{
public:
	
	/** @brief Number of equations */
	int fNumEq;
	/** @brief Number of Threads */
	int fNumthreads;
	/** @brief Polynomial order */
	int fPolyOrder;
	/** @brief Number of Substructures */
	int fNumSub;
	/** @brief Number of coarse equations */
	int fNumEqCoarse;
	/** @brief Number of iterations */
	int fniter;
	/** @brief Number of elements in fMultiply */
	int fnMultiply;
	/** @brief Number of elements in fPreCond */
	int fnPreCond;
	/** @brief Number of elements in the mesh */
	int fNumberofElements;
	/** @brief Time for Substructuring Mesh */
	REAL ft0sub;
	/** @brief Time for Computing the system of equations for each substructure */
	REAL ft1comput;
	/** @brief Time for Convert Graph */
	REAL ft2congraph;
	/** @brief Time for AnalyseGraph */
	REAL ft3analysegraph;
	/** @brief Total Time for Identifying Corner Nodes = Convert Graph + AnalyseGraph */
	REAL ft4identcorner;
	/** @brief Time for ThreadDohrmanAssembly */
	REAL ft5dohrassembly;
	/** @brief Time to decompose the inner nodes matrix */
	REAL ft55decompmatriznosinternos;
	/** @brief Total Time for Iterations */
	REAL ft6iter;
	
	/** @brief Time to Multiply for each iteration */
	TPZStack<REAL> fMultiply;
	/** @brief Time to PreCond for each iteration */
	TPZStack<REAL> fPreCond;
	
	/** @brief Initialize the attributes */
	TPZTimeTemp();
	
	/** @brief Append to the file a line with the informations */
	void PrintLine(std::ostream &out);
	
	/** @brief Read one line of the file */
	void ReadLine(std::istream &in);
	
	/** @brief Print the header */
	static void PrintHeader(std::ostream &out);
	
	/** @brief If the file does not exists, returns "true" */
	static bool NeedsHeader(std::string &FileName);
	
private:
	
};

/**
 * @ingroup util
 */
/// External variable to TPZTimeTemp (to take time)
extern TPZTimeTemp tempo;


#endif
