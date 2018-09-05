/**
 * @file
 * @brief Contains the TPZPoint class which defines the topology of a point. 
 */

#ifndef PZTOPOLOGYTPZPOINT_H
#define PZTOPOLOGYTPZPOINT_H

#include "pzreal.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pztrnsform.h"
#include "pzeltype.h"
#include "pzstack.h"

class TPZIntPoints;
class TPZInt1Point;
class TPZGraphEl1dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a point. \ref topology "Topology"
	 * It has a one side (the same element).
	 */
	class TPZPoint : public TPZSavable {
	public:
		
		/** @brief Enumerate for topological characteristics */
		enum {NCornerNodes = 1, NSides = 1, Dimension = 0, NFaces = 0};

                virtual int ClassId() const;
                void Read(TPZStream& buf, void* context);
                void Write(TPZStream& buf, int withclassid) const;

                
		/** @brief Default constructor */
        TPZPoint() : TPZRegisterClassId(&TPZPoint::ClassId) {
		}
		/** @brief Default destructor */
		virtual ~TPZPoint() {
		}

		/** @name About sides of the topological element
		 * @{ */
		
		/** @brief Returns the dimension of the side */
		static int SideDimension(int side) {
			return 0;
		}

		/** @brief Get all sides with lower dimension on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides) {
		}
		/** @brief Get all sides with lower dimension but equal to DimTarget on side */
		static void LowerDimensionSides(int side,TPZStack<int> &smallsides, int targetdim) {
		}
		
		/**
		 * @brief Returns all sides whose closure contains side
		 * @param side Smaller dimension side
		 * @param high Vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high) {
		}
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side) {
			return 1;
		}
		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node) {
			return 0;
		}

		/** @brief Returns number of connects of the element  ??? */
		static int NumSides() { return 1; }
		/** @brief Returns the number of connects for a set dimension   // Jorge ??? */
		static int NumSides(int dimension) { return 0; };	

		/** @brief Returns the number of nodes (not connectivities) associated with a side   // Jorge - sides or nodes??? */
		static int NContainedSides(int side) {
			return 1;
		}
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c) {
			return 0;
		}

		/** @} */

		/** @name About points at the parametric spaces
		 * @{ */

		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center) {
		}

		/** @brief Verifies if the parametric point pt is in the element parametric domain */
		static bool IsInParametricDomain(TPZVec<REAL> &pt, REAL tol = 1e-6){
			return true;
		}
        
        /** @brief Generates a random point in the master domain */
        static void RandomPoint(TPZVec<REAL> &pt)
        {
            
        }
        
        template<class T>
        static bool MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide);
        
        static void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);

		/** @} */

		/** @name About type of the topological element
		 * @{ */

		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static MElementType Type();
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side) ;

		/** @} */
		
		/** @name About Transformations
		 * @{ */

		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom side where the point resides
		 * @param sideto side whose closure contains sidefrom
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto) {
			TPZTransform<> result(0,0);
			return result;
		}
		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 * @see the class TPZTransform
		 */
		static TPZTransform<> TransformSideToElement(int side) {
			TPZTransform<> result(0,0);
			return result;
		}

		static TPZTransform<> TransformElementToSide(int side){
			TPZTransform<> t(0,0);
			return t;
		}
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int64_t> &id);
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs
		 * of the corner nodes
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(int side, TPZVec<int64_t> &id);
		
		/** @} */

		/** @name Methods related over numeric integration
		 * @{ */
		
		/**
		 * @brief Create an integration rule over side
		 * @param side Side to create integration rule
		 * @param order Order of the integration rule to be created
		 */
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZInt1Point IntruleType;
		
		/* @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather)
	{
		permgather[0] = 0;
		return;
	}
		
		/** @brief Volume of the master element (measure of the element) */
		static REAL RefElVolume() {
			return 0.;
		}
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao)
        {
            sides[0] = 0;
            dir[0] = 0;
            bilinearounao[0] = 0;
        }

        static void GetSideDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors)
        {
            sides[0] = 0;
            dir[0] = 0;
            bilinearounao[0] = 0;
            sidevectors[0] = 0;
        }

        /// Compute the directions of the HDiv vectors
        static void ComputeDirections(TPZFMatrix<REAL> &gradx, REAL detjac, TPZFMatrix<REAL> &directions)
        {
        }
        

        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();

	};
	
}

#endif
