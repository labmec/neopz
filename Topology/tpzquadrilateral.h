/**
 * @file
 * @brief Contains the TPZQuadrilateral class which defines the topology of a quadrilateral element. 
 */
// C++ Interface: tpzquadrilateral
//
// Description: 
//
//
// Author: Philippe R. B. Devloo <phil@fec.unicamp.br>, (C) 2007
//
// Copyright: See COPYING file that comes with this distribution
//
//
#ifndef PZTOPOLOGYTPZQUADRILATERAL_H
#define PZTOPOLOGYTPZQUADRILATERAL_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"

#ifdef _AUTODIFF
#include "fadType.h"
#endif

class TPZIntPoints;
class TPZIntQuad;
class TPZGraphElQ2dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a quadrilateral element. \ref topology "Topology"
	 */
	class TPZQuadrilateral {
	public:
		/// Topological data
		enum {NSides = 9, NCornerNodes = 4, Dimension = 2};
		/// Simple constructor
		TPZQuadrilateral();
		/// Destructor
		virtual ~TPZQuadrilateral();
		
		/// get all sides with lower dimension on side
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/// get all sides with lower dimension but equal to DimTarget on side
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		/**
		 * @brief return the number of nodes (not connectivities) associated with a side
		 */
		static int NSideNodes(int side);
		/**
		 * @brief returns the local node number of the node "node" along side "side"
		 */
		static int SideNodeLocId(int side, int node);
		/**
		 * @brief returns the barycentric coordinates in the master element space of the original element
		 */
		static void CenterPoint(int side, TPZVec<REAL> &center);
		/** @brief volume of the master element */
		static REAL RefElVolume(){return 4.0;}
		
		/**
		 * @brief Create an integration rule 
		 * @param order order of the integration rule to be created
		 * @param side side to create integration rule
		 */
		static TPZIntPoints * CreateSideIntegrationRule(int side, int order);
		
		/// To numerical integration rule
		typedef TPZIntQuad IntruleType;
		/// To graphical element type
		typedef TPZGraphElQ2dd GraphElType;
		
		/// return the type of the element as specified in file pzeltype.h
		static MElementType Type();// { return EQuadrilateral;}
		
		/**
		 * @brief return the type of the element as specified in file pzeltype.h
		 */
		static MElementType Type(int side);
		
		
		/**
		 * @brief Number of connects of the element (9)
		 * @return number of connects of the element
		 */
		static int NumSides()
		{
			return NSides;
		}
		
		/**
		 * @brief returns the transformation which takes a point from the side sidefrom ot
		 * the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 */
		static TPZTransform SideToSideTransform(int sidefrom, int sideto);
		
		/**
		 * @brief returns the dimension of the side
		 */
		static int SideDimension(int side);
		
		/**
		 * @brief return the number of nodes (not connectivities) associated with a side
		 */
		static int NContainedSides(int side);
		/**
		 * @brief return the number of connects for a set dimension
		 */
		static int NumSides(int dimension);
		/**
		 * @brief returns the local connect number of the connect "c" along side "side"
		 */
		static int ContainedSideLocId(int side, int c);
		/**
		 * @brief return the connect associate to side side is a particular method for hdiv space
		 **/
		//static int ContainedSideLocId(int side);
		
		/**
		 * @brief Returns the transformation which transform a point from the interior of the element to the side
		 * @param side side to which the point will be tranformed (0<=side<=8)
		 * @return TPZTransform object
		 * @see the class TPZTransform
		 */
		static TPZTransform TransformElementToSide(int side);
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform object
		 * @see the class TPZTransform
		 */
		static TPZTransform TransformSideToElement(int side);
		
		/** @brief Verifies if the parametric point pt is in the element parametric domain
		 */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);
		/// function pointer which determines the type of computational element
		/**
		 * function pointer which determines what type of computational element will be created
		 * Method which identifies the transformation based on the IDs
		 * of the corner nodes
		 * @param id indexes of the corner nodes
		 * @return index of the transformation of the point corresponding to the topology
		 */
		//static TPZCompEl *(*fp)(TPZGeoEl *el,TPZCompMesh &mesh,int &index);
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
