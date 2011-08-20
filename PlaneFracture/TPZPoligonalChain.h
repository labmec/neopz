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

using namespace std;

#include <set>
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzvec.h"

/** 
 * @brief Implements a poligonal chain 
 * @author Cesar Lucci
 * @since 29/09/2010
 */
class TPZPoligonalChain
{
	public:
	
	TPZPoligonalChain();
	~TPZPoligonalChain();
	
			/** @brief Data for Dot of the monotone chain */ 
			struct MCDot
			{
				public:
				
				MCDot()
				{
					fCoord.Resize(2);
					SetId(-1);
					SetCoordX(0.);
					SetCoordY(0.);
				}
				
				void SetId(int Id)
				{
					fId = Id;
				}
				void SetCoordX(double x)
				{
					fCoord[0] = x;
				}
				void SetCoordY(double y)
				{
					fCoord[0] = y;
				}
				
				int GetId()
				{
					return fId;
				}
				double GetCoordX()
				{
					return fCoord[0];
				}
				double GetCoordY()
				{
					return fCoord[1];
				}
				
				int fId;
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
				std::map<double, MCDot> GetDotMap()
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


	static int Signal(double val)
	{
		int sig = -1;
		if(val >= 0.)
		{
			sig = +1;
		}		
		return sig;
	}
	
	void SplitInMonotoneChains(TPZVec<REAL> &PoligonalChain, TPZFMatrix &FromR3toR2);
	bool ThereIsIntersection(MCDot &p0, MCDot &p1, MCDot &q0, MCDot &q1, MCDot &intersec);
	
	//	Continuar daqui!!!
	//	> fazer o metodo MCI (do mathematica)
	
	int NMCs()
	{
		return fMCVec.size();
	}
	
	MC GetMC(int pos)
	{
		#ifdef DEBUG
		if(pos >= NMCs())
		{
			std::cout << "Invalid Monotone Chain on " << __PRETTY_FUNCTION__ << std::endl;
			DebugStop();
		}
		#endif
		
		return fMCVec[pos];
	}
	
	std::vector<MC> fMCVec;
};

#endif