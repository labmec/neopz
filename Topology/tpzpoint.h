/**
 * @file
 * @brief Contains the TPZPoint class which defines the topology of a point. 
 */
// C++ Interface: tpzpoint
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZPOINT_H
#define PZTOPOLOGYTPZPOINT_H

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "tpzpoint.h"
#include "pzeltype.h"
#include "tpzint1point.h"

class TPZIntPoints;
class TPZInt1Point;
class TPZGraphEl1dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a point. \ref topology "Topology"
	 */
	class TPZPoint {
	public:
		
		enum {NCornerNodes = 1, NSides = 1, Dimension = 0};
		
		typedef TPZInt1Point IntruleType;
		
		TPZPoint();
		
		virtual ~TPZPoint();
		
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides) {
			//    smallsides.Push(0);
		}
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int targetdim) {
			//    smallsides.Push(0);
		}
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high) {
		}
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NSideNodes(int side) {
			return 1;
		}
		/**
		 * @brief Returns the local node number of the node "node" along side "side"
		 */
		static int SideNodeLocId(int side, int node) {
			return 0;
		}
		/**
		 * @brief Returns the barycentric coordinates in the master element space of the original element
		 */
		static void CenterPoint(int side, TPZVec<REAL> &center) {
		}
		/** @brief volume of the master element*/
		static REAL RefElVolume(){
			return 0.;
		}
		/**
		 * @brief Returns the dimension of the side
		 */
		static int SideDimension(int side) {
			return 0;
		}
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom ot
		 * the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 */
		static TPZTransform SideToSideTransform(int sidefrom, int sideto) {
			TPZTransform result(0,0);
			return result;
		}
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform object
		 * @see the class TPZTransform
		 */
		static TPZTransform TransformSideToElement(int side) {
			TPZTransform result(0,0);
			return result;
		}
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain
		 */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6){
			return true;
		}
		
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static MElementType Type() ;//{ return EPoint;}
		
		static std::string StrType() { return "Point";}
		
		/**
		 * @brief Returns the type of the element as specified in file pzeltype.h
		 */
		static MElementType Type(int side) ;
		
		static std::string StrType(int side) { return "Point";}
		
		
		static int NumSides() {return 1;}
		/**
		 * @brief Returns the number of nodes (not connectivities) associated with a side
		 */
		static int NContainedSides(int side) {
			return 1;
		}
		/**
		 * @brief Returns the number of connects for a set dimension
		 */
		static int NumSides(int dimension){return 0;};	
		/**
		 * @brief Returns the local connect number of the connect "c" along side "side"
		 */
		static int ContainedSideLocId(int side, int c) {
			return 0;
		}
		
		static TPZTransform TransformElementToSide(int side){
			TPZTransform t(0,0);
			return t;
		}
		
		/**
		 * @brief Function pointer which determines what type of computational element will be created
		 */
		static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
		/**
		 * @brief Method which identifies the transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int> &id);
		
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
