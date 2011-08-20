/**
 * \file
 * @brief Contains implementations of the TPZPoligonalChain methods.
 */
/*
 *  TPZPoligonalChain.cpp
 *  Crack
 *
 *  Created by Cesar Lucci on 29/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "TPZPoligonalChain.h"

TPZPoligonalChain::TPZPoligonalChain()
{	
	fMCVec.resize(0);
}

TPZPoligonalChain::~TPZPoligonalChain()
{
}

void TPZPoligonalChain::SplitInMonotoneChains(TPZVec<REAL> &PoligonalChain, TPZFMatrix &FromR3toR2)
{
	#ifdef DEBUG
	// a poligonal chain is defined by points in R3,
	// so PoligonalChain vector must have size multiple of 3 {x0, y0, z0, x1, y1, z1, x2, y2, z2...}
	if(PoligonalChain.NElements()%3 != 0)
	{
		std::cout << "Invalid Poligonal Chain on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	int nPtos = PoligonalChain.NElements()/3;
	
	int MCsign = +1, dXsign;
	double dX, coordX0 = 0., coordX1, coordY1, X, Y, Z;
	TPZFMatrix coordR3(3,1), coordR2(2,1);
	
	MC actualMC;
	MCDot actualDot;
	
	for(int p = 0; p < nPtos; p++)
	{
		X = PoligonalChain[3*p];
		Y = PoligonalChain[3*p+1];
		Z = PoligonalChain[3*p+2];
		
		coordR3(0,0) = X;
		coordR3(1,0) = Y;
		coordR3(2,0) = Z;
		
		FromR3toR2.Multiply(coordR3, coordR2);
		coordX1 = coordR2(0,0);
		coordY1 = coordR2(1,0);
		
		dX = coordX1 - coordX0;
		dXsign = Signal(dX);
		if(dXsign != MCsign)
		{
			if(actualMC.fDotMap.size() > 0)
			{
				fMCVec.push_back(actualMC);
				actualMC.Reset();
				actualMC.fDotMap[coordX0] = actualDot;	
			}
			MCsign = Signal(dX);		
		}
		
		actualDot.SetId(p);		
		actualDot.SetCoordX(coordX1);
		actualDot.SetCoordY(coordY1);
		
		actualMC.GetDotMap()[coordX1] = actualDot;
		
		coordX0 = coordX1;
	}
}

bool TPZPoligonalChain::ThereIsIntersection(MCDot &p0, MCDot &p1, MCDot &q0, MCDot &q1, MCDot &intersec)
{
	bool thereIs = false;
	
	// λ e μ sao solucoes do sistema: (1-λ).p0+(λ).p1 == (1-μ).q0+(μ).q1
	// se (0 ≤ λ ≤ 1) e (0 ≤ μ ≤ 1), as semiretas p0→p1 e q0→q1 interceptam-se!
	
	double p0x = p0.GetCoordX();
	double p0y = p0.GetCoordY();
	double p1x = p1.GetCoordX();
	double p1y = p1.GetCoordY();
	
	double q0x = q0.GetCoordX();
	double q0y = q0.GetCoordY();
	double q1x = q1.GetCoordX();
	double q1y = q1.GetCoordY();
	
	double numerador, denominador, lambda, mu;
	
	numerador = -(q0y*q1x) + p0y*(-q0x + q1x) + p0x*(q0y - q1y) + q0x*q1y;
	denominador = p1y*(q0x - q1x) + p0y*(-q0x + q1x) + (p0x - p1x)*(q0y - q1y);
	if(fabs(denominador) < 1.E-15)
	{
		return thereIs;
	}
	lambda = numerador/denominador;
	
	numerador = p0y*(p1x - q0x) + p1y*q0x - p1x*q0y + p0x*(-p1y + q0y);
	denominador = p1y*(q0x - q1x) + p0y*(-q0x + q1x) + (p0x - p1x)*(q0y - q1y);
	if(fabs(denominador) < 1.E-15)
	{
		return thereIs;
	}
	mu = numerador/denominador;
	
	if(lambda >= 0. && lambda <= 1. && mu >= 0. && mu <= 1.)
	{
		thereIs = true;
		
		double intersecX = (1.-lambda)*p0x + lambda*p1x;
		double intersecY = (1.-lambda)*p0y + lambda*p1y;
		
		intersec.SetCoordX(intersecX);
		intersec.SetCoordY(intersecY);
	}
	
	return thereIs;
}