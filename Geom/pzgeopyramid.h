/**
 * @file
 * @brief Contains the TPZGeoPyramid class which implements the geometry of pyramid element.
 */

#ifndef TPZGEOTETRAPIRAMIDH
#define TPZGEOTETRAPIRAMIDH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "tpzpyramid.h"

template<class TVar>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of pyramid element. \ref geometry "Geometry"
	 */
	class TPZGeoPyramid  : public TPZNodeRep<5, pztopology::TPZPyramid>
	{
	public:
		/** @brief Number of corner nodes */
		enum {NNodes = 5};
		
		/** @brief Constructor with list of nodes */
		TPZGeoPyramid(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoPyramid() : TPZNodeRep<NNodes, pztopology::TPZPyramid>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoPyramid(const TPZGeoPyramid &cp,
					  std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPyramid(const TPZGeoPyramid &cp) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPyramid(const TPZGeoPyramid &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Pyramid";}
		
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
		static  void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		/** @brief Computes the geometric location*/
		static  void X(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZVec<REAL> &result);
		
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
										  TPZVec<long>& nodeindexes,
										  int matid,
										  long& index);
	};
	
};

#endif 
