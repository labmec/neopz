/**
 * @file
 * @brief Contains the TPZLine class which defines the topology of a line element.
 */

#ifndef PZTOPOLOGYTPZLINE_H
#define PZTOPOLOGYTPZLINE_H

#include "pzfmatrix.h"
#include "pzstack.h"
#include "pztrnsform.h"
#include "pzquad.h"
#include "pzeltype.h"
#include "pzaxestools.h"

#include "fadType.h"

#include "TPZTopologyUtils.h"

class TPZIntPoints;
class TPZInt1d;
class TPZGraphEl1dd;

class TPZCompEl;
class TPZGeoEl;
class TPZCompMesh;

/// Groups all classes defining the structure of the master element
namespace pztopology {
	/**
	 * @ingroup topology
	 * @author Philippe R. B. Devloo
	 * @brief Defines the topology of a line element. \ref topology "Topology"
	 * Sides 0 and 1 are vertices, side 2 is the line. 
	 */
	class TPZLine : public TPZSavable {
	public:
    friend void pztopology::GetPermutation<TPZLine>(const int permute, TPZVec<int> &permutation);
		/** @brief Topological characteristics */
    static constexpr int64_t NSides = 3;
    static constexpr int64_t NCornerNodes = 2;
    static constexpr int64_t Dimension = 1;
    static constexpr int64_t NFacets = 2;
    static constexpr int64_t NPermutations = 2;
      
		
    int ClassId() const override;
    void Read(TPZStream &buf, void *context) override;
    void Write(TPZStream &buf, int withclassid) const override;

		/** @brief Default constructor */
        TPZLine() : TPZRegisterClassId(&TPZLine::ClassId){
		}
		
		/** @brief Default destructor */
		virtual ~TPZLine() {
		}

        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(const TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);

		/** @name About sides of the topological element
		 * @{ */

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

		/** @brief Returns the number of connects for a set dimension */
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
		static constexpr MElementType Type() {return EOned;}
		/** @brief Returns the type of the element side as specified in file pzeltype.h */
		static MElementType Type(int side);

		/** @} */
		
		/** @name About Transformations
		 * @{ */

		/**
		 * @brief Returns the transformation which takes a point from the side sidefrom to the side sideto
		 * @param sidefrom Side where the point resides
		 * @param sideto Side whose closure contains sidefrom
		 */
		static TPZTransform<> SideToSideTransform(int sidefrom, int sideto);

		/**
		 * @brief Returns the transformation which transform a point from the side to the interior of the element
		 * @param side Side from which the point will be tranformed (0<=side<=2)
		 * @return TPZTransform<> object
		 */
		static TPZTransform<> TransformSideToElement(int side);

		/**
		 * @brief Returns the transformation which transform a point from the interior of the element to the side
		 * @param side Side to which the point will be tranformed (0<=side<=2)
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
		 * @brief Method which identifies the transformation of a side based on the IDs
		 * of the corner nodes
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
		typedef TPZInt1d IntruleType;
		/** @brief Typedef to graphical element type */
		typedef TPZGraphEl1dd GraphElType;
		
		/** @} */
		
		/**
		 * @brief Identifies the permutation of the nodes needed to make neighbouring elements compatible 
		 * in terms of order of shape functions
		 * @param side Side for which the permutation is needed
		 * @param id Ids of the corner nodes of the elements
		 * @param permgather Permutation vector in a gather order
		 */
		static void GetSideHDivPermutation(int transformationid, TPZVec<int> &permgather);

		/** @brief Volume of the master element (measure of the element) */
		static constexpr REAL RefElVolume() { return 2.0; }
        
        /* Given side and gradx the method returns directions needed for Hdiv space */
        static void ComputeDirections(int side, TPZFMatrix<REAL> &gradx, TPZFMatrix<REAL> &directions, TPZVec<int> &sidevectors);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao);
        static void GetSideHDivDirections(TPZVec<int> &sides, TPZVec<int> &dir, TPZVec<int> &bilinearounao, TPZVec<int> &sidevectors);
        
        /// Compute the directions of the HDiv vectors
        template <class TVar>
        static void ComputeHDivDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions)
        {
            TVar detjac = TPZAxesTools<TVar>::ComputeDetjac(gradx);
            
            for (int i=0; i<3; i++) {
                directions(i,0) = -gradx(i,0)/detjac;
                directions(i,1) = gradx(i,0)/detjac;
                directions(i,2) = gradx(i,0)/detjac;
            }
        }

        template <class TVar>
        static void ComputeHCurlDirections(TPZFMatrix<TVar> &gradx, TPZFMatrix<TVar> &directions, const TPZVec<int> &transformationIds){
            const REAL sign = transformationIds[0] == 0 ? 1 : -1;
            directions(0,0) = 0.5 * sign;
            directions(0,1) = 0.5 * sign;
            directions(0,2) = 0.5 * sign;
        }
        
        /// Compute the directions of the HDiv vectors for constant divergent
        // template <class TVar>
        static void ComputeConstantHDiv(TPZVec<REAL> &point, TPZFMatrix<REAL> &vecDiv, TPZVec<REAL> &div){
            DebugStop();
        };
        static void ComputeConstantHCurl(TPZVec<REAL> &point, TPZFMatrix<REAL> &vecDiv, TPZFMatrix<REAL> &curl, const TPZVec<int> &transformationIds){
            DebugStop();
        };
        
        /**
         * Returns the number of bilinear sides to this shape. Needed to compute the number shapefunctions( NConnectShapeF )
         */
        static int NBilinearSides();
	protected:
        /** @brief Valid permutations between nodes*/
        static constexpr int fPermutations [2][3] ={
          {0,1,2},
          {1,0,2}
        };
		
	};
	
}

#endif
