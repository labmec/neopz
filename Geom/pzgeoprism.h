// $Id: pzgeoprism.h,v 1.15 2011-05-11 01:38:41 phil Exp $
/**
 * @file
 * @brief Contains the TPZGeoPrism class which implements the geometry of a prism element.
 */

#ifndef TPZGEOPRISMH
#define TPZGEOPRISMH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzprism.h"

template <class TVar>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

#include <string>


namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a prism element. \ref geometry "Geometry"
	 */
	class TPZGeoPrism : public TPZNodeRep<6, pztopology::TPZPrism>  
	{
	public:
		
		enum {NNodes = 6};
		/** @brief Constructor with list of nodes */
		TPZGeoPrism(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZPrism>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoPrism() : TPZNodeRep<NNodes, pztopology::TPZPrism>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoPrism(const TPZGeoPrism &cp,
					std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPrism(const TPZGeoPrism &cp) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPrism(const TPZGeoPrism &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZPrism>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Prism";} 
		
		/** @brief Implementation of two-dimensional bilinear interpolation*/
		static  void Shape(TPZVec<REAL> &x,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);

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
        
		/** @brief Computes the jacobian*/
		static  void Jacobian(TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		/** @brief Computes the geometric location*/
		static  void X(TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);
		
		/** @brief Returns the projection of a given point from "NSide - 1" side to "side". */
		static bool MapToSide(int side, TPZVec<REAL> &InternalPar, TPZVec<REAL> &SidePar, TPZFMatrix<REAL> &JacToSide);
		
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
		 */
		/**
		 * This points are CornerNodes, when projected in the opposing side
		 */
		static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);
		
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid, int& index);
	};
	
};
#endif 
