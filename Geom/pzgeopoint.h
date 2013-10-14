/**
 * @file
 * @brief Contains the TPZGeoPoint class which implements the geometry of a point element.
 */

#ifndef TPZGEOPOINTH
#define TPZGEOPOINTH

#include "pznoderep.h"
#include "pzvec.h"
#include "pzeltype.h"
#include "pzintel.h"
#include "tpzpoint.h"

template<class TVar> 
class TPZFMatrix;
class TPZGeoEl;
class TPZIntPoints;
class TPZInt1Point;
class TPZGraphEl1dd;
class TPZGeoMesh;

#include <string>


/**
 * @brief Groups all classes which model the geometry
 * @see TPZIntelGen
 */
/** 
 * Objects of this class implement the mapping between the master element
 * and deformed element
 * These classes are used as template arguments of @see TPZGeoElement and
 */
namespace pzgeom {
	
	/** 
	 * @ingroup geometry
	 * @brief Implements the geometry of a point element. \ref geometry "Geometry"
	 */
	class TPZGeoPoint : public TPZNodeRep<1, pztopology::TPZPoint> {
		
	public:
		enum {NNodes = 1};
        
		/** @brief Auxiliar structure to accellerate computations */
		struct TMem {
		};

        typedef pztopology::TPZPoint Top;
		
		/** @brief Constructor with list of nodes */
		TPZGeoPoint(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZPoint>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoPoint() : TPZNodeRep<NNodes, pztopology::TPZPoint>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoPoint(const TPZGeoPoint &cp,
					std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZPoint>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPoint(const TPZGeoPoint &cp) : TPZNodeRep<NNodes, pztopology::TPZPoint>(cp)
		{
		}
		
		/** @brief Copy constructor with given mesh */
		TPZGeoPoint(const TPZGeoPoint &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZPoint>(cp)
		{
		}
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Point";}
        
        void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
		void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }

		static void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi);
		
		static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian) {
			jacobian.Redim(nodes.Rows(),0);
		}
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
	public:
		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid, long& index);
	};
	
};

#endif

