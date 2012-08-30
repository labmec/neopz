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
}

TPZPoligonalChain::~TPZPoligonalChain()
{
}

void TPZPoligonalChain::SplitInMonotoneChains(TPZVec<REAL> &PoligonalChain, std::vector<MC> &MCVec)
{
    MCVec.resize(0);
    
	#ifdef DEBUG
	// a poligonal chain is defined by points in R3,
	// so PoligonalChain vector must have size multiple of 3 {x0, y0, z0, x1, y1, z1, x2, y2, z2...}
	if(PoligonalChain.NElements()%2 != 0)
	{
		std::cout << "Invalid Poligonal Chain on " << __PRETTY_FUNCTION__ << std::endl;
		DebugStop();
	}
	#endif
	
	int nPtos = PoligonalChain.NElements()/2;
	
	int MCsign = +1, dXsign;
	double dX, coordX0 = 0., coordX1, coordZ1;
	
	MC actualMC;
	MCDot actualDot;
	
	for(int p = 0; p < nPtos; p++)
	{
		coordX1 = PoligonalChain[2*p];
		coordZ1 = PoligonalChain[2*p+1];
		
		dX = coordX1 - coordX0;
		dXsign = Signal(dX);
		if(dXsign != MCsign)
		{
			if(actualMC.fDotMap.size() > 0)
			{
				MCVec.push_back(actualMC);
				actualMC.Reset();
				actualMC.fDotMap[coordX0] = actualDot;	
			}
			MCsign = dXsign;		
		}
		
		actualDot.SetCoordX(coordX1);
		actualDot.SetCoordZ(coordZ1);
		
		actualMC.GetDotMap()[coordX1] = actualDot;
		
		coordX0 = coordX1;
	}
}

void TPZPoligonalChain::JustSplitAndPrint(TPZVec<REAL> &PoligonalChain, std::ostream & out)
{
    std::vector<MC> MCVec(0);
    SplitInMonotoneChains(PoligonalChain, MCVec);
    
    int nMC = MCVec.size();
    for(int mc = 0; mc < nMC; mc++)
    {
        out << "Monotone Chain #" << mc << ":\n";
        MC actMC = MCVec[mc];
        
        int npts = actMC.NDots();
        for(int p = 0; p < npts; p++)
        {
            out << "x = " << actMC.GetDot(p).GetCoordX() << " | z = " << actMC.GetDot(p).GetCoordZ() << std::endl;
        }
        
        out << "\n------------------------------\n";
    }
    out.flush();
}

bool TPZPoligonalChain::ThereIsIntersection(MCDot &p0, MCDot &p1, MCDot &q0, MCDot &q1, MCDot &intersec)
{
	bool thereIs = false;
	
	// λ e μ sao solucoes do sistema: (1-λ).p0+(λ).p1 == (1-μ).q0+(μ).q1
	// se (0 ≤ λ ≤ 1) e (0 ≤ μ ≤ 1), as semiretas p0→p1 e q0→q1 interceptam-se!
	
	double p0x = p0.GetCoordX();
	double p0z = p0.GetCoordZ();
	double p1x = p1.GetCoordX();
	double p1z = p1.GetCoordZ();
	
	double q0x = q0.GetCoordX();
	double q0z = q0.GetCoordZ();
	double q1x = q1.GetCoordX();
	double q1z = q1.GetCoordZ();
	
	double numerador, denominador, lambda, mu;
	
	numerador = -(q0z*q1x) + p0z*(-q0x + q1x) + p0x*(q0z - q1z) + q0x*q1z;
	denominador = p1z*(q0x - q1x) + p0z*(-q0x + q1x) + (p0x - p1x)*(q0z - q1z);
	if(fabs(denominador) < 1.E-15)
	{
		return thereIs;
	}
	lambda = numerador/denominador;
	
	numerador = p0z*(p1x - q0x) + p1z*q0x - p1x*q0z + p0x*(-p1z + q0z);
	denominador = p1z*(q0x - q1x) + p0z*(-q0x + q1x) + (p0x - p1x)*(q0z - q1z);
	if(fabs(denominador) < 1.E-15)
	{
		return thereIs;
	}
	mu = numerador/denominador;
	
	if(lambda >= 0. && lambda <= 1. && mu >= 0. && mu <= 1.)
	{
		thereIs = true;
		
		double intersecX = (1.-lambda)*p0x + lambda*p1x;
		double intersecZ = (1.-lambda)*p0z + lambda*p1z;
		
		intersec.SetCoordX(intersecX);
		intersec.SetCoordZ(intersecZ);
	}
	
	return thereIs;
}