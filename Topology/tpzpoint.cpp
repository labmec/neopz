/**
 * @file
 * @brief Contains the implementation of the TPZPoint methods. 
 */
//
// C++ Implementation: tpzpoint
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#include "tpzpoint.h"
#include "pzquad.h"
#include "tpzint1point.h"
#include "pzeltype.h"

#include "pzcreateapproxspace.h"

namespace pztopology {
	
	TPZCompEl *(*TPZPoint::fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index) = CreatePointEl;
	
	
	TPZPoint::TPZPoint()
	{
	}
	
	
	TPZPoint::~TPZPoint()
	{
	}
	
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
	
	/**
	 * Method which identifies the transformation based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */
	int TPZPoint::GetTransformId(TPZVec<int> &id)
	{
		return 0;
	}
	
	/**
	 * Method which identifies the transformation of a side based on the IDs
	 * of the corner nodes
	 * @param id indexes of the corner nodes
	 * @return index of the transformation of the point corresponding to the topology
	 */	
	int TPZPoint::GetTransformId(int side, TPZVec<int> &id)
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
	void TPZPoint::GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather)
	{
		permgather[0] = 0;
		return;
	}
	
}
