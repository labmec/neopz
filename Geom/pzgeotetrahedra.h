/**
 * @file
 * @brief Contains the TPZGeoTetrahedra class which implements the geometry of a tetrahedral element.
 */

#ifndef TPZGEOTETRAHEDRAH
#define TPZGEOTETRAHEDRAH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpztetrahedron.h"

#include <string>

template<class TVar>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/** 
	 * @ingroup geometry
	 * @brief Implements the geometry of a tetrahedral element. \ref geometry "Geometry"
	 */
	class TPZGeoTetrahedra  : public TPZNodeRep<4,pztopology::TPZTetrahedron>
	{
	public:
		
		enum {NNodes = 4};
		/** @brief Constructor with list of nodes */
		TPZGeoTetrahedra(TPZVec<int> &nodeindexes) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoTetrahedra() : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp,
						 std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp, TPZGeoMesh &) : TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Tetra";} 
		
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
		
		/** @brief Returns the projection of a given point from \f$ NSide - 1 \f$ side to \f$ side \f$. */
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
		/** This points are CornerNodes, when projected in the opposing side */
		static void FixSingularity(int side, TPZVec<REAL>& OriginalPoint, TPZVec<REAL>& ChangedPoint);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int>& nodeindexes,
										  int matid,
										  int& index);
	};
	
};

#endif 
