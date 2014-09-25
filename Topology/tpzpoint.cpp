/**
 * @file
 * @brief Contains the implementation of the TPZPoint methods. 
 */

#include "tpzpoint.h"
#include "pzquad.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"

namespace pztopology {
	
	TPZIntPoints *TPZPoint::CreateSideIntegrationRule(int side, int order)
	{
		return new IntruleType(order);
	}
	
	
	MElementType TPZPoint::Type()
	{
		return EPoint;
	}
	
	MElementType TPZPoint::Type(int side)
	{
		switch(side) {
			case 0:
				return EPoint;
			default:
				return ENoType;
		}
	}
    
    int TPZPoint::NBilinearSides()
    {return 0;}
    
    bool TPZPoint::MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide) {
		SidePar.Resize(0); JacToSide.Resize(0,0);
		return true;
	}
    
    void TPZPoint::ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord)
    {
        //It doesnt make any sense ask node coordinate of 0D element
        DebugStop();
    }
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZPoint::GetTransformId(TPZVec<long> &id)
	{
		return 0;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZPoint::GetTransformId(int side, TPZVec<long> &id)
	{
		return 0;
	}
	
	/**
	 * Identifies the permutation of the nodes needed to make neighbouring elements compatible 
	 * in terms of order of shape functions
	 * @param side : side for which the permutation is needed
	 * @param id : ids of the corner nodes of the elements
	 * @param permgather : permutation vector in a gather order
	 */
//	void TPZPoint::GetSideHDivPermutation(int side, TPZVec<long> &id, TPZVec<int> &permgather)
//	{
//		permgather[0] = 0;
//		return;
//	}
    void TPZPoint::ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors)
    {
        directions.Redim(3, 0);
        sidevectors.Resize(0);
    }
	
}
