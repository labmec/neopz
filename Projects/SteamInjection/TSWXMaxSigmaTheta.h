//$Id: TSWXMaxSigmaTheta.h,v 1.5 2009-07-03 19:43:07 caju Exp $

/*
 *  TSWXMaxSigmaTheta.h
 *  TSWXMaxSigmaTheta
 *
 *  Created by CAJU on 6/18/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 *  Esta classe calcula o SigmaTheta m·ximo em um poÁo vertical, reservatÛrio circular.
 *  Neste c·lculo È utilizada uma malha bidimensional axi-simÈtrica para modelar o reservatÛrio
 *  e o SigmaTheta m·ximo È localizado no ponto (rw,0), ou seja, na interface entre poÁo-reservatÛrio no fundo do poÁo.
 */

#ifndef MaxSigmaThetaH
#define MaxSigmaThetaH

#include "pzreal.h"
#include <iostream>
#include <cstdlib>
#include <math.h>

#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeopoint.h"
#include "pzcmesh.h"
#include "pzmaterial.h"
#include "pzelasAXImat.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"
#include "pzmaterialdata.h"
#include "pzgengrid.h"

using namespace std;
using namespace pzgeom;

//-- DO NOT TOUCH -------------------------------------------------------- FROM HERE
//elemento
const int matElId = 1;

//contorno
const int matLeftDOWNid = -1;
const int matLeftUPid = -2;
const int matRightDOWNid = -3;
const int matRightUPid = -4;
const int matTOPid = -5;
const int matBOTTOMid = -6;
const int matMiddleUPid = -7;
const int matMiddleDOWNid = -8;

//bc's types
const int dirichlet = 0;
const int newmann = 1;
const int mista = 2;
const int hidrostatica = 3;
//------------------------------------------------------------------------------- TO HERE

class TSWXMaxSigmaTheta
{
public:
	
	TSWXMaxSigmaTheta();

	/**
	 * @param distr : value of distributed load (might be injection pressure? - ask Phil, Im using unit pressure)
 	 * @param rw: well radius
	 * @param h : heigh of distributed load = heigh of reservoir (just porous matrix)
	 * @param E - Young Modulus
	 * @param nu - Poison ratio
	 */
	TSWXMaxSigmaTheta(double rw, double h, double E, double nu);
	~TSWXMaxSigmaTheta();
	void SetData(double rw, double h, double E, double nu);
	
	//b is the position of distributed load
	void ComputeMaxSigmaTheta(double t, double b, double DistrRightDown);
	
	void PrintMathematica(ostream &out = std::cout);

	map<double,double> & GetSolution()
	{
		return fb_sol;
	}

	double getE() { return fE; }
	double getNu() { return fnu; }
	
private:

	TPZGeoMesh * GetGeoMesh()
	{
		return fGmesh;
	}
	TPZCompMesh * GetCompMesh()
	{
				return fCmesh;
	}
	
	//b is the position of distributed load	
	void CreateGeoMesh(double b);

	//DistrRightDown is the value of pressure in the vapor front
	void CreateCompMesh(double DistrRightDown);
	
	TPZGeoMesh * fGmesh;
	TPZCompMesh * fCmesh;
	
	//Directional refinement intensity
	int fDirectNdiv;
	
	//Uniform refinement intensity
	int fUnifNDiv;
	
	//Maximum quantity of Mesh subdivisions (maximum refinement rate)
	int fDivMAX;
	
	//interpolation order
	int fOrder;
	
	//Geometry
	double frw;	//well radius
	double fB;		//2D reservoir width
	double fH;		//2D reservoir heigh
	double fh;		//heigh of distributed load
	
	//Materials
	double fE;		//Young modulus
	double fnu;	//Poisson ratio
	
	//Applied Loads
	double fDistrLeftDown;		//Carregamento distribuido no contorno aa esquerda, entre 0 e h
	double fDistrLeftUp;			//Carregamento distribuido no contorno aa Esquerda, entre h e H
	double fDistrMiddleDown;	//Carregamento distribuido numa linha vertical interna aa malha, entre 0 e h
	double fDistrMiddleUp;		//Carregamento distribuido numa linha vertical interna aa malha, entre h e H
	double fDistrRightUp;			//Carregamento distribuido no contorno aa Direita, entre h e H
	double fDistrSurface;			//Carregamento distribuido no interior da malha
	
	//Solution, i.e., desired output (b vs Maximum SigmaTheta)
	map<double,double> fb_sol;
};
#endif
