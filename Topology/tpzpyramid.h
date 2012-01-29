/**
 * @file
 * @brief Contains the TPZPyramid class which defines the topology of a Pyramid element. 
 */
// C++ Interface: tpzpyramid
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZPYRAMID_H
#define PZTOPOLOGYTPZPYRAMID_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZIntPyram3D;
class TPZGraphElPyramidMapped;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a Pyramid element. \ref topology "Topology"
	 */
	class TPZPyramid {
	public:
		
		enum {NSides = 19, NCornerNodes = 5, Dimension = 3};
		
		TPZPyramid();
		
		virtual ~TPZPyramid();
		
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NSideNodes(int side);
		/**
		 * @brief Returns the local node number of the node "node" along side "side"
		 */
		static int SideNodeLocId(int side, int node);
		/**
		 * @brief Returns the barycentric coordinates in the master element space of the original element
		 */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		/** @brief Volume of the master element*/
		static REAL RefElVolume(){return (4./3.);}
		/**
		 * @brief Returns the dimension of the side
		 */
		static int SideDimension(int side);
		/**
		 * @brief Returns the number of connects for a set dimension
		 */
		static int NumSides(int dimension);
		
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom ot
		 * the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 */
		static TPZTransform SideToSideTransform(int sidefrom, int sideto);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side side to which the point will be tranformed (0<=side<=20)
		 * @return TPZTransform object
		 * @see the class TPZTransform
		 */
		static TPZTransform TransformElementToSide(int side);
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side side from which the point will be tranformed (0<=side<=20)
		 * @return TPZTransform object
		 * @see the class TPZTransform
		 */
		static TPZTransform TransformSideToElement(int side);
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain
		 */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);
		
		/**
		 * @brief Create an integration rule 
		 * @param order order of the integration rule to be created
		 * @param side side to create integration rule
		 */
		static TPZIntPoints * CreateSideIntegrationRule(int side, int order);
		
		typedef TPZIntPyram3D IntruleType;
		typedef TPZGraphElPyramidMapped GraphElType;
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static MElementType Type();// { return EPiramide;}
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static MElementType Type(int side);
		
		/**
		 * @brief Number of connects of the element (21)
		 * @return number of connects of the element
		 */
		static int NumSides();
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NContainedSides(int side);
		/**
		 * @brief Returns the local connect number of the connect "c" along side "side"
		 */
		static int ContainedSideLocId(int side, int c);
		
		/**
		 * @brief Function pointer which determines what type of computational element will be created
		 */
		//static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
		
	protected:
		/** 
		 * @brief Data structure which defines the pyramid transformations and topology
		 */
		static int FaceNodes[5][4];
		/** 
		 * @brief Data structure which defines the pyramid transformations and topology
		 */
		static int SideNodes[8][2];
		/** 
		 * @brief Data structure which defines the pyramid transformations and topology
		 */
		static int ShapeFaceId[5][4];
		
	public:
		/**
		 * @brief Method which identifies the transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int> &id)
        {
            DebugStop();
            return -1;
        }
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs
		 * of the corner nodes
		 * @param side index of side
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(int side, TPZVec<int> &id);
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side : side for which the permutation is needed
		 * @param id : ids of the corner nodes of the elements
		 * @param permgather : permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather);
		
	};
	
}

#endif
