/**
 * @file
 * @brief Contains the TPZGeoCube class which implements the geometry of hexahedra element.
 */
// $Id: TPZGeoCube.h,v 1.12 2011-05-11 01:38:40 phil Exp $

#ifndef TPZGEOCUBEH
#define TPZGEOCUBEH


#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpzcube.h"

#include <string>

template<class TVar>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of hexahedra element. \ref geometry "Geometry"
	 */
	class TPZGeoCube : public TPZNodeRep<8, pztopology::TPZCube> {
		
	public:
		enum {NNodes = 8};
		
		/** @brief Constructor with list of nodes */
		TPZGeoCube(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZCube>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoCube() : TPZNodeRep<NNodes, pztopology::TPZCube>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoCube(const TPZGeoCube &cp,
				   std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoCube(const TPZGeoCube &cp) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoCube(const TPZGeoCube &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZCube>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Hexa";} 
		
		/* brief compute the coordinate of a point given in parameter space */
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        /* @brief compute the jacobian of the map between the master element and deformed element */
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }
        
		static void X(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		/** @brief Returns the projection of a given point from \f$ NSide - 1 \f$ side to \f$ side \f$. */
		static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void Jacobian(TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
							 TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
		
	};
	
};

#endif
