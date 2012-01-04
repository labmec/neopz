/**
 * @file
 * @brief Contains the TPZPoligonalChain class which defines a poligonal chain.
 */
#ifndef TPZPOLIGONALCHAINH
#define TPZPOLIGONALCHAINH

/*
 *  TPZPoligonalChain.h
 *  Crack
 *
 *  Created by Cesar Lucci on 29/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <set>
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzvec.h"

#include <ostream>

/** 
 * @brief Implements a poligonal chain 
 * @author Cesar Lucci
 * @since 29/09/2010
 */

/** @brief Data for Dot of the monotone chain */ 
struct MCDot
{
    public:
    
    MCDot()
    {
        fCoord.Resize(2);
        SetCoordX(0.);
        SetCoordZ(0.);
    }
    
    void SetCoordX(double x)
    {
        fCoord[0] = x;
    }
    
    void SetCoordZ(double z)
    {
        fCoord[1] = z;
    }
    
    double GetCoordX()
    {
        return fCoord[0];
    }
    
    double GetCoordZ()
    {
        return fCoord[1];
    }
    
    TPZVec<REAL> fCoord;
};

/** @brief Monotone chain */
struct MC
{
    public:
    
    MC()
    {
        Reset();
    }
    
    void Reset()
    {
        fDotMap.clear();
    }
    
    int NDots()
    {
        return fDotMap.size();
    }
    
    std::map<double, MCDot> & GetDotMap()
    {
        return fDotMap;
    }
    
    MCDot GetDot(int pos)
    {
        #ifdef DEBUG
        if(pos >= NDots())
        {
            std::cout << "Invalid Dot on " << __PRETTY_FUNCTION__ << std::endl;
            DebugStop();
        }
        #endif
        
        MCDot DotFound;
        
        std::map<double, MCDot>::iterator it;
        int p = 0;
        for(it = fDotMap.begin(); it != fDotMap.end(); it++)
        {
            if(p == pos)
            {
                DotFound = it->second;
            }
            p++;
        }
        return DotFound;
    }
    
    std::map<double, MCDot> fDotMap; // indexed by x value
};

//---------------------------------------------------------------

class TPZPoligonalChain
{
	public:
	
	TPZPoligonalChain();
	~TPZPoligonalChain();
	
	static int Signal(double val)
	{
		int sig = -1;
		if(val >= 0.)
		{
			sig = +1;
		}		
		return sig;
	}
	
	static void SplitInMonotoneChains(TPZVec<REAL> &PoligonalChain, std::vector<MC> &MCVec);
    static void JustSplitAndPrint(TPZVec<REAL> &PoligonalChain, std::ostream & out = std::cout);
	static bool ThereIsIntersection(MCDot &p0, MCDot &p1, MCDot &q0, MCDot &q1, MCDot &intersec);
	
	//	Continuar daqui!!!
	//	> fazer o metodo MCI (do mathematica)
};

#endif