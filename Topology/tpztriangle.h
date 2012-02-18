/**
 * @file
 * @brief Contains the TPZTriangle class which defines the topology of a triangle. 
 */

#ifndef PZTOPOLOGYTPZTRIANGLE_H
#define PZTOPOLOGYTPZTRIANGLE_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZGraphElTd;
class TPZIntTriang;
class TPZGraphElT2dMapped;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @brief Defines the topology of a triangle element. \ref topology "Topology"
	 * Sides 0 to 2 are vertices, sides 3 to 5 are lines, side 6 is the triangle. 
	 * @author Philippe R. B. Devloo
	 */
	class TPZTriangle {
	public:
		
		/** @brief Enumerate for topological characteristics */
		enum {NSides = 7, NCornerNodes= 3, Dimension = 2};
		
		/** @brief Default constructor */
		TPZTriangle(){
		};
		
		/** @brief Default destructor */
		virtual ~TPZTriangle(){
		};
		
		/** @name About sides of the topological element
		 * @{ */
		
		/** @brief Returns the dimension of the side */
		static int SideDimension(int side);

		/** @brief Get all sides with lower dimension on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides);
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int DimTarget);
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side smaller dimension side
		 * @param high vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);

		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side);
		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);
		
		/** @brief Returns the number of connects of the element (7) */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NContainedSides(int side);
		/** @brief Returns the local connect along side "side" especial for hdivspace */
		static int ContainedSideLocId(int side);
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);

		/** @} */
		
		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);

		/** @brief Verifies if the parametric point pt is in the element parametric domain */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6);

		/** @} */
		
		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type();// { return ETriangle;}
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type(int side);
		
		/** @} */
		
		/** @name About Transformations
		 * @{ */
		
		/** 
		 * @brief Returns the transformation which takes a point from the side sidefrom of the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 * @see the class TPZTransform
		 */
		static TPZTransform SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform object
		 */
		static TPZTransform TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which transform a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=8)
		 * @return TPZTransform object
		 */
		static TPZTransform TransformElementToSide(int side);
		
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int> &id);
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @param side Side of the element in which the transformation is
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(int side, TPZVec<int> &id);
		
		/** @} */
				
		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints * CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZIntTriang IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElT2dMapped GraphElType;

		/** @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather);
		
		/** @brief Volume of the master element (measure) */
		static REAL RefElVolume() { return 0.5L; }

	};
	
}

#endif
