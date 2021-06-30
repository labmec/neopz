/**
 * @file
 * @brief Contains the TPZPyramid class which defines the topology of a pyramid element. 
 */

#ifndef PZTOPOLOGYTPZPYRAMID_H
#define PZTOPOLOGYTPZPYRAMID_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "pzaxestools.h"
#include "TPZTopologyUtils.h"

class TPZIntPoints;
class TPZIntPyram3D;
class TPZGraphElPyramidMapped;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {

	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a Pyramid element. \ref topology "Topology"
	 * Sides 0 to 4 are vertices, sides 5 to 12 are lines, side 13 are quadrilateral (pyramid base),
	 * sides 14 to 17 are triangles and side 18 is the pyramid.
	 */
	class TPZPyramid : public TPZSavable{
	public:
    friend void pztopology::GetPermutation<TPZPyramid>(const int permute, TPZVec<int> &permutation);
		/** @brief Topological characteristics */
		static constexpr int64_t NSides = 19;
    static constexpr int64_t NCornerNodes = 5;
    static constexpr int64_t Dimension = 3;
    static constexpr int64_t NFacets = 5;
    static constexpr int64_t NPermutations = 8;
      
		
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;

		/** @brief Default constructor */
        TPZPyramid() : TPZRegisterClassId(&TPZPyramid::ClassId) {
		}
		
		/** @brief Default destructor */
		virtual ~TPZPyramid() {
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
		
		/** @brief Number of connects of the element (21) */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NContainedSides(int side);
		/** @brief Returns the local connect number of the connect "c" along side "side" */
		static int ContainedSideLocId(int side, int c);


        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        /**
         * This method calculates the influence (a.k.a. the blend function) of the side side regarding an
         * interior point qsi. It is used by the TPZGeoBlend class.
         * @param side the index of the side
         * @param xi coordinates of the interior point
         * @param blendFactor influence (0 <= blendFactor <= 1)
         * * @param corrFactorDxi derivative of the blendFactor in respect to xi
         */
        template<class T>
        static void BlendFactorForSide(const int &side, const TPZVec<T> &xi, T &blendFactor,
                                      TPZVec<T> &corrFactorDxi);
		/** @} */

		/** @name About points at the parametric spaces
		 * @{ */
				
		/** @brief Returns the barycentric coordinates in the master element space of the original element */
		static void CenterPoint(int side, TPZVec<REAL> &center);

        /** @brief Verifies if the parametric point pt is in the element parametric domain */
        static bool IsInParametricDomain(const TPZVec<REAL> &pt, REAL tol = pztopology::gTolerance);

        /** @brief Verifies if the parametric point pt is in the element parametric domain (FAD version)*/
        static bool IsInParametricDomain(const TPZVec<Fad<REAL>> &pt, REAL tol = pztopology::gTolerance){
            TPZVec<REAL> xi(pt.size());
            for(int i = 0; i < pt.size(); i++) xi[i]= pt[i].val();
            return IsInParametricDomain(xi,tol);
        }
        /** @brief Generates a random point in the master domain */
        static void RandomPoint(TPZVec<REAL> &pt);

        /**
         * This method will check if the projection to a certain side (MapToSide method) is regular,
         * i.e., if the interior point in the parametric domain is not too close to the projection's singularity.
         * @param side the index of the side upon which the interior point will be projected upon
         * @param xiInterior coordinates of the interior point
         * @return true if the interior point is far from the singularity
         */
        template<class T>
        static bool CheckProjectionForSingularity(const int &side, const TPZVec<T> &xiInterior);

        template<class T>
        static void MapToSide(int side, TPZVec<T> &InternalPar, TPZVec<T> &SidePar, TPZFMatrix<T> &JacToSide);
        
        static void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
		
		/** @} */

		/** @name About type of the topological element
		 * @{ */
		
		/** @brief Returns the type of the element as specified in file pzeltype.h */
		static constexpr MElementType Type() {return EPiramide;}
		
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
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=20)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=20)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(TPZVec<int64_t> &id)
        {
            DebugStop();
            return -1;
        }
		
		/**
		 * @brief Method which identifies the transformation of a side based on the IDs of the corner nodes
		 * @param side Index of side
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */	
		static int GetTransformId(const int side, const TPZVec<int64_t> &id);

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
		typedef TPZIntPyram3D IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElPyramidMapped GraphElType;

		/** @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather);
		
		/** @brief Volume of the master element*/
		static constexpr REAL RefElVolume() {return (4.L/3.L); }
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        template <class TVar>
        static void ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions);
        
        /// Adjust the directions associated with the tip of the pyramid, considering that one of the faces is constrained
        template <class TVar>
        static void AdjustTopDirections(int ConstrainedFace,TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions);
        

        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();
        
        /**
         * @brief Computes the corner shape functions of the element
         * @param pt (input) point where the shape function is computed
         * @param phi (output) value of the (5) shape functions
         * @param dphi (output) value of the derivatives of the (5) shape functions holding the derivatives in a column
         */
        static void CornerShape(const TPZVec<REAL> &pt, TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

	protected:
		/** @name Data structure which defines the pyramid transformations */
		/** @{ */
    /** @brief Nodes over quadrilateral sides (2d - faces). */
    static constexpr int FaceNodes[5][4]  = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
        
    /** @brief Nodes over lines sides (1d) */
    static constexpr int SideNodes[8][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,4},{2,4},{3,4} };

		/** @brief Ids of the shape face */
    static constexpr int ShapeFaceId[5][4] = { {0,1,2,3},{0,1,4,-1},{1,2,4,-1},{3,2,4,-1},{0,3,4,-1} };
	
		/** @brief Valid permutations between nodes*/
    static constexpr int fPermutations[8][19] = {
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18},/*00*/
      {0,3,2,1,4,8,7,6,5,9,12,11,10,13,17,16,15,14,18},/*01*/
      {1,0,3,2,4,5,8,7,6,10,9,12,11,13,14,17,16,15,18},/*02*/
      {1,2,3,0,4,6,7,8,5,10,11,12,9,13,15,16,17,14,18},/*03*/
      {2,1,0,3,4,6,5,8,7,11,10,9,12,13,15,14,17,16,18},/*04*/
      {2,3,0,1,4,7,8,5,6,11,12,9,10,13,16,17,14,15,18},/*05*/
      {3,0,1,2,4,8,5,6,7,12,9,10,11,13,17,14,15,16,18},/*06*/
      {3,2,1,0,4,7,6,5,8,12,11,10,9,13,16,15,14,17,18} /*07*/
    };
                /** @} */
		
	};
	
}

#endif
