// -*- c++ -*-

// $ Id: $

// TPZGeoQuad.h: interface for the TPZGeoQuad class.
//
//////////////////////////////////////////////////////////////////////

#ifndef TPZGEOQUADH
#define TPZGEOQUADH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzquadrilateral.h"

#include <string>


class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/// @ingroup geometry
	/// @brief Implements the geometry of a quadrilateral element
	class TPZGeoQuad  : public TPZNodeRep<4, pztopology::TPZQuadrilateral>
	{
	public:
		
		enum {NNodes = 4};
		enum {NVectors = 18};
		/**
		 * @brief Constructor with list of nodes
		 */
		TPZGeoQuad(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(nodeindexes)
		{
		}
		
		/**
		 * @brief Empty constructor
		 */
		TPZGeoQuad() : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>()
		{
		}
		
		/**
		 * @brief Constructor with node map
		 */
		TPZGeoQuad(const TPZGeoQuad &cp,
				   std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp,gl2lcNdMap)
		{
		}
		
		/**
		 * @brief Copy constructor
		 */
		TPZGeoQuad(const TPZGeoQuad &cp) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/**
		 * @brief Copy constructor
		 */
		TPZGeoQuad(const TPZGeoQuad &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZQuadrilateral>(cp)
		{
		}
		
		/**
		 * @brief Returns the type name of the element
		 */
		static std::string TypeName() { return "Quad";} 
		
		/** @brief Implementation of two-dimensional bilinear interpolation*/
		static  void Shape(TPZVec<REAL> &x,TPZFMatrix &phi,TPZFMatrix &dphi);
		
		/** @brief Implementation of normal vector to Hdiv space*/
		/** 
		 construct the normal vector for element Hdiv
		 */
		static void VecHdiv(TPZFMatrix & coord,TPZFMatrix &fNormalVec,TPZVec<int> & fVectorSide);
		/** @brief Computes the vecorial product of the two vectors*/ 
		static void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);
		/** @brief Computes normal vector to plane determinated by three points */
		static void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);
		
		/** @brief Computes the jacobian*/
		static  void Jacobian(TPZFMatrix & coord, TPZVec<REAL>& par, TPZFMatrix &jacobian, TPZFMatrix &axes,REAL &detjac,TPZFMatrix &jacinv);
		
		/** @brief Computes the geometric location*/
		static  void X(TPZFMatrix & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);
		
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
