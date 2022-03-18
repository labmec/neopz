/**
 * @file
 * @brief Contains the TPZCube class which defines the topology of the hexahedron element.
 */

#ifndef PZTOPOLOGYTPZCUBE_H
#define PZTOPOLOGYTPZCUBE_H


#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzeltype.h"
#include "pzaxestools.h"

#include "fadType.h"

class TPZIntPoints;
class TPZIntCube3D;
class TPZGraphElQ3dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;
#include "TPZTopologyUtils.h"
/// Groups all classes defining the structure of the master element
namespace pztopology {
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of the hexahedron element. \ref topology "Topology"
	 * Sides 0 to 7 are vertices, sides 8 to 19 are lines, 20 to 25 are quadrilaterals 
	 * and side 26 is the hexahedra (cube).
	 */
	class TPZCube : public TPZSavable {
	public:
    friend void pztopology::GetPermutation<TPZCube>(const int permute, TPZVec<int> &permutation);
		/** @brief Topological characteristics */
    static constexpr int64_t NSides = 27;
    static constexpr int64_t NCornerNodes = 8;
    static constexpr int64_t Dimension = 3;
    static constexpr int64_t NFacets = 6;
    static constexpr int64_t NPermutations = 48;
      
		
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;
                
		/** @brief Default constructor */
        TPZCube() : TPZRegisterClassId(&TPZCube::ClassId) {
		}
		
		/** @brief Default destructor */
		virtual ~TPZCube() {
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
		 * @param side Smaller dimension side
		 * @param high Vector which will contain all sides whose closure contain sidefrom
		 */
		static void HigherDimensionSides(int side, TPZStack<int> &high);
		
		/** @brief Returns the number of nodes (not connectivities) associated with a side */
		static int NSideNodes(int side);
		/** @brief Returns the local node number of the node "node" along side "side" */
		static int SideNodeLocId(int side, int node);

		/** @brief Returns number of connects of the element (27)  ??? */
		static int NumSides();
		/** @brief Returns the number of connects for a set dimension */
		static int NumSides(int dimension);
		
		/** @brief Returns the number of connectivities associated with a side */
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
		static constexpr MElementType Type() { return ECube;}
		
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
		 * @param side Side from which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);
		/**
		 * @brief Returns the transformation which projects a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=26)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformElementToSide(int side);
		
		/**
		 * @brief Method which identifies the transformation based on the IDs of the corner nodes
		 * @param id Indexes of the corner nodes
		 * @return Index of the transformation of the point corresponding to the topology
		 */
		static int GetTransformId(const TPZVec<int64_t> &id);
		
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
		static TPZIntPoints *CreateSideIntegrationRule(int side, int order);
		
		/** @brief Typedef to numerical integration rule */
		typedef TPZIntCube3D IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphElQ3dd GraphElType;

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
		static constexpr REAL RefElVolume(){return 8.0;}
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        template <class TVar>
        static void ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions);
        
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);

        /// Compute the directions of the HDiv vectors for constant divergent
        // template <class TVar>
        static void ComputeConstantHDiv(TPZVec<REAL> &point, TPZFMatrix<REAL> &vecDiv, TPZVec<REAL> &div);
        static void ComputeConstantHCurl(TPZVec<REAL> &point, TPZFMatrix<REAL> &vecDiv, TPZFMatrix<REAL> &curl, const TPZVec<int> &transformationIds);
        static int GetSideOrient(const int &face);
        
        /** Compute the directions of the HCurl vectors.
         * These vectors are combined with H1 shape functions to create the HCurl shape functions.
         * They *must be* computed in the following order:
         * - \f$v^{e,a}\f$: vector associated with edge \f$e\f$. It is normal to the edge \f$\hat{e}\f$ adjacent to \f$e\f$e by the vertex \f$a\f$a.
         * - \f$v^{e,T}\f$: vector associated with edge \f$e\f$. It is tangent to the edge \f$\hat{e}\f$.
         * - \f$v^{F,e}\f$: vector associated with face \f$F\f$. It is normal to the face \f$\hat{F}\f$ adjacent to \f$F\f$e by the edge \f$e\f$a.
         * - \f$v^{F,T}\f$: two orthornormal vectors associated with face \f$F\f$ and tangent to it.
         * - \f$v^{F,\perp}\f$: outward normal vector associated with face \f$F\f$ (3D only)
         * - \f$v^{K}\f$: set of orthonormal vectors associated with the volume of the element itself (3D only. In 2D \f$v^{F,T}\f$ does its job)
         * The side ordering should be respected. In the definition of the \f$v^{e,a}\f$ and the \f$v^{F,e}\f$ vectors, the subsides are ordered as the return of LowerDimensionSides.
         * @tparam TVar REAL or Fad<REAL>
         * @param gradx the gradient of the element mapping. if computing in normal element, gradx is the identity matrix.
         * @param directions computed directions
         * @param transformationIds transformation Ids associated with each side of dim > 0
         */
        template <class TVar>
        static void ComputeHCurlDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions, const TPZVec<int> &transformationIds);
        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();
	
	protected:
		/** @name Data structure which defines the hexahedral transformations */
		/** @{ */
		/** @brief For each face was enumerated the pontoal sides (vertices) */
    static constexpr int FaceNodes[6][4]  = { {0,1,2,3},{0,1,5,4},{1,2,6,5},{3,2,6,7},{0,3,7,4},{4,5,6,7} };
	
	/** @brief For each edge was enumerated the pontoal sides (vertices) */
	static constexpr int SideNodes[12][2]  = { {0,1},{1,2},{2,3},{3,0},{0,4},{1,5},
		{2,6},{3,7},{4,5},{5,6},{6,7},{7,4} };
	
	/** @brief For each face was enumerated the vertice sides on its main diagonal */
	static constexpr int ShapeFaceId[6][2] = { {0,2},{0,5},{1,6},{3,6},{0,7},{4,6} };

    /** @brief Valid permutations between nodes*/
    static constexpr int fPermutations[48][27] = {
      {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26},/*000*/
      {0,1,5,4,3,2,6,7,8,13,16,12,11,9,17,19,10,14,18,15,21,20,22,25,24,23,26},/*001*/
      {0,3,2,1,4,7,6,5,11,10,9,8,12,15,14,13,19,18,17,16,20,24,23,22,21,25,26},/*002*/
      {0,3,7,4,1,2,6,5,11,15,19,12,8,10,18,16,9,14,17,13,24,20,23,25,21,22,26},/*003*/
      {0,4,5,1,3,7,6,2,12,16,13,8,11,19,17,9,15,18,14,10,21,24,25,22,20,23,26},/*004*/
      {0,4,7,3,1,5,6,2,12,19,15,11,8,16,18,10,13,17,14,9,24,21,25,23,20,22,26},/*005*/
      {1,0,3,2,5,4,7,6,8,11,10,9,13,12,15,14,16,19,18,17,20,21,24,23,22,25,26},/*006*/
      {1,0,4,5,2,3,7,6,8,12,16,13,9,11,19,17,10,15,18,14,21,20,24,25,22,23,26},/*007*/
      {1,2,3,0,5,6,7,4,9,10,11,8,13,14,15,12,17,18,19,16,20,22,23,24,21,25,26},/*008*/
      {1,2,6,5,0,3,7,4,9,14,17,13,8,10,18,16,11,15,19,12,22,20,23,25,21,24,26},/*009*/
      {1,5,4,0,2,6,7,3,13,16,12,8,9,17,19,11,14,18,15,10,21,22,25,24,20,23,26},/*010*/
      {1,5,6,2,0,4,7,3,13,17,14,9,8,16,18,10,12,19,15,11,22,21,25,23,20,24,26},/*011*/
      {2,1,0,3,6,5,4,7,9,8,11,10,14,13,12,15,17,16,19,18,20,22,21,24,23,25,26},/*012*/
      {2,1,5,6,3,0,4,7,9,13,17,14,10,8,16,18,11,12,19,15,22,20,21,25,23,24,26},/*013*/
      {2,3,0,1,6,7,4,5,10,11,8,9,14,15,12,13,18,19,16,17,20,23,24,21,22,25,26},/*014*/
      {2,3,7,6,1,0,4,5,10,15,18,14,9,11,19,17,8,12,16,13,23,20,24,25,22,21,26},/*015*/
      {2,6,5,1,3,7,4,0,14,17,13,9,10,18,16,8,15,19,12,11,22,23,25,21,20,24,26},/*016*/
      {2,6,7,3,1,5,4,0,14,18,15,10,9,17,19,11,13,16,12,8,23,22,25,24,20,21,26},/*017*/
      {3,0,1,2,7,4,5,6,11,8,9,10,15,12,13,14,19,16,17,18,20,24,21,22,23,25,26},/*018*/
      {3,0,4,7,2,1,5,6,11,12,19,15,10,8,16,18,9,13,17,14,24,20,21,25,23,22,26},/*019*/
      {3,2,1,0,7,6,5,4,10,9,8,11,15,14,13,12,18,17,16,19,20,23,22,21,24,25,26},/*020*/
      {3,2,6,7,0,1,5,4,10,14,18,15,11,9,17,19,8,13,16,12,23,20,22,25,24,21,26},/*021*/
      {3,7,4,0,2,6,5,1,15,19,12,11,10,18,16,8,14,17,13,9,24,23,25,21,20,22,26},/*022*/
      {3,7,6,2,0,4,5,1,15,18,14,10,11,19,17,9,12,16,13,8,23,24,25,22,20,21,26},/*023*/
      {4,0,1,5,7,3,2,6,12,8,13,16,19,11,9,17,15,10,14,18,21,24,20,22,25,23,26},/*024*/
      {4,0,3,7,5,1,2,6,12,11,15,19,16,8,10,18,13,9,14,17,24,21,20,23,25,22,26},/*025*/
      {4,5,1,0,7,6,2,3,16,13,8,12,19,17,9,11,18,14,10,15,21,25,22,20,24,23,26},/*026*/
      {4,5,6,7,0,1,2,3,16,17,18,19,12,13,14,15,8,9,10,11,25,21,22,23,24,20,26},/*027*/
      {4,7,3,0,5,6,2,1,19,15,11,12,16,18,10,8,17,14,9,13,24,25,23,20,21,22,26},/*028*/
      {4,7,6,5,0,3,2,1,19,18,17,16,12,15,14,13,11,10,9,8,25,24,23,22,21,20,26},/*029*/
      {5,1,0,4,6,2,3,7,13,8,12,16,17,9,11,19,14,10,15,18,21,22,20,24,25,23,26},/*030*/
      {5,1,2,6,4,0,3,7,13,9,14,17,16,8,10,18,12,11,15,19,22,21,20,23,25,24,26},/*031*/
      {5,4,0,1,6,7,3,2,16,12,8,13,17,19,11,9,18,15,10,14,21,25,24,20,22,23,26},/*032*/
      {5,4,7,6,1,0,3,2,16,19,18,17,13,12,15,14,8,11,10,9,25,21,24,23,22,20,26},/*033*/
      {5,6,2,1,4,7,3,0,17,14,9,13,16,18,10,8,19,15,11,12,22,25,23,20,21,24,26},/*034*/
      {5,6,7,4,1,2,3,0,17,18,19,16,13,14,15,12,9,10,11,8,25,22,23,24,21,20,26},/*035*/
      {6,2,1,5,7,3,0,4,14,9,13,17,18,10,8,16,15,11,12,19,22,23,20,21,25,24,26},/*036*/
      {6,2,3,7,5,1,0,4,14,10,15,18,17,9,11,19,13,8,12,16,23,22,20,24,25,21,26},/*037*/
      {6,5,1,2,7,4,0,3,17,13,9,14,18,16,8,10,19,12,11,15,22,25,21,20,23,24,26},/*038*/
      {6,5,4,7,2,1,0,3,17,16,19,18,14,13,12,15,9,8,11,10,25,22,21,24,23,20,26},/*039*/
      {6,7,3,2,5,4,0,1,18,15,10,14,17,19,11,9,16,12,8,13,23,25,24,20,22,21,26},/*040*/
      {6,7,4,5,2,3,0,1,18,19,16,17,14,15,12,13,10,11,8,9,25,23,24,21,22,20,26},/*041*/
      {7,3,0,4,6,2,1,5,15,11,12,19,18,10,8,16,14,9,13,17,24,23,20,21,25,22,26},/*042*/
      {7,3,2,6,4,0,1,5,15,10,14,18,19,11,9,17,12,8,13,16,23,24,20,22,25,21,26},/*043*/
      {7,4,0,3,6,5,1,2,19,12,11,15,18,16,8,10,17,13,9,14,24,25,21,20,23,22,26},/*044*/
      {7,4,5,6,3,0,1,2,19,16,17,18,15,12,13,14,11,8,9,10,25,24,21,22,23,20,26},/*045*/
      {7,6,2,3,4,5,1,0,18,14,10,15,19,17,9,11,16,13,8,12,23,25,22,20,24,21,26},/*046*/
      {7,6,5,4,3,2,1,0,18,17,16,19,15,14,13,12,10,9,8,11,25,23,22,21,24,20,26} /*047*/
    };
		/** @} */
        
		
	};
	
}

#endif
