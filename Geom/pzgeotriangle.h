// -*- c++ -*-
// $ Id: $
// TPZGeoQuad.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOTRIANGLEH
#define TPZGEOTRIANGLEH

#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpztriangle.h"

#include <string>
#include <map>

class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;


namespace pzgeom {
	
	/// @ingroup geometry
	/// @brief Implements the geometry of a triangle element
	class TPZGeoTriangle : public TPZNodeRep<3, pztopology::TPZTriangle> 
	{
	public:
		
		enum {NNodes = 3};
		
		
		/**
		 * @brief Constructor with list of nodes
		 */
		TPZGeoTriangle(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes)
		{
		}
		
		/**
		 * @brief Empty constructor
		 */
		TPZGeoTriangle() : TPZNodeRep<NNodes,pztopology::TPZTriangle>()
		{
		}
		
		/**
		 * @brief Constructor with node map
		 */
		TPZGeoTriangle(const TPZGeoTriangle &cp,
					   std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap)
		{
		}
		
		/**
		 * @brief Copy constructor
		 */
		TPZGeoTriangle(const TPZGeoTriangle &cp) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}
		
		/**
		 * @brief Copy constructor
		 */
		TPZGeoTriangle(const TPZGeoTriangle &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}
		
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "Triangle";} 
		
		/** @brief Implementation of two-dimensional bilinear interpolation*/
		static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		/** @brief Computes the jacobian*/
		static  void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);
		
		/** @brief Computes the geometric location*/
		static  void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);
		/** @brief Implementation of Hdiv space*/
		static	void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);
		
		static	void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);
		
		static void VecHdiv(TPZFMatrix & coord, TPZFMatrix & fNormalVec,TPZVec<int> &sidevector);
		
		
		
		/**
		 * @brief Returns the projection of a given point from "NSide - 1" side to "side".
		 */
		static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix &JacToSide);
		
		/**
		 * @brief Method which creates a geometric boundary condition 
		 * element based on the current geometric element, 
		 * a side and a boundary condition number
		 */
		static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
	protected:
		/**
		 * @brief This method apply an infinitesimal displacement in some points
		 * to fix singularity problems when using MapToSide() method!
		 * This points are CornerNodes, when projected in the opposing side
		 */
		static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);
		
	public:
		/**
		 * @brief Creates a geometric element according to the type of the father element
		 */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
		
	};
	
};

#endif 
