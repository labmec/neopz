/**
 * @file
 * @brief Contains the TPZPrism class which defines the topology of a Prism. 
 */

#ifndef PZTOPOLOGYTPZPRISM_H
#define PZTOPOLOGYTPZPRISM_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"

class TPZIntPoints;
class TPZIntPrism3D;
class TPZGraphElPrismMapped;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a Prism. \ref topology "Topology"
	 * Sides 0 to 7 are vertices, sides 7 to 14 are lines, sides 15 and 19 are triangles, 
	 * sides 16 to 18 are quadrilaterals and side 20 is the prism.
	 */
	class TPZPrism {
	public:
					
		/** @brief Enumerate for topological characteristics */
		enum {NSides = 21, NCornerNodes = 6, Dimension = 3};
		
		/** @brief Default constructor */
		TPZPrism() {
		}
		
		/** @brief Default destructor */
		virtual ~TPZPrism() {
		}
		
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
		
		/** @brief Returns the number of connects of the element (21) */
		static int NumSides();
		
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NContainedSides(int side);
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);

		/** @} */
		
		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);

		/** @brief Verifies if the parametric point pt is in the element parametric domain */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6L);
		
		/** @} */
		
		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type();
		
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side);

		/** @} */
		
		/** @name About Transformations
		 * @{ */
		
		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 * @see the class TPZTransform
		 */
		static TPZTransform SideToSideTransform(int sidefrom, int sideto);
		
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=20)
		 * @return TPZTransform object
		 */
		static TPZTransform TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=20)
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
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
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
		typedef TPZIntPrism3D IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElPrismMapped GraphElType;

		/** @} */

		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int side, TPZVec<int> &id, TPZVec<int> &permgather);

		/** @brief Volume of the master element*/
		static REAL RefElVolume() { return 1.0L; }
		
	protected:
		/** @name Data structure which defines the prism transformations and topology */
		/** @{ */
		
		/** @brief Nodes over quadrilateral sides (2d - faces). */
		static int FaceNodes[5][4];

		/** @brief Nodes over lines sides (1d) */
		static int SideNodes[9][2];

		/** @brief Ids of the shape face */
		static int ShapeFaceId[5][4];

		/** @} */
		
	};
	
}

#endif
