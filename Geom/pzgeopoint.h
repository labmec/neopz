/**
 * @file
 * @brief Contains the TPZGeoPoint class which implements the geometry of a point element or 0-D element.
 */

#ifndef TPZGEOPOINTH
#define TPZGEOPOINTH

#include <map>            // for map, operator==
#include <sstream>        // for basic_stringbuf<>::int_type, basic_stringbu...
#include "pzeltype.h"     // for MElementType
#include "pzmatrix.h"     // for TPZFMatrix, TPZMatrix, TPZFNMatrix
#include "pznoderep.h"    // for TPZNodeRep
#include "pznoderep.h.h"  // for TPZNodeRep::TPZNodeRep<N, Topology>
#include "pzreal.h"       // for REAL
#include "tpzpoint.h"     // for TPZPoint
class TPZGeoEl;
class TPZGeoMesh;
template <class T> class TPZVec;


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
        
        static bool IsLinearMapping(int side)
        {
            return true;
        }
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Point";}
        
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
		
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            gradx.Zero();
        }
        
        template<class T>
		static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result);
        
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
		
		static void Shape(TPZVec<REAL> &pt,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi)
        {
            phi(0,0) = 1.;
        }
		
        template<class T>
        static void TShape(TPZVec<T> &pt,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi)
        {
            phi(0,0) = (T)1.;
        }
        
		static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
		
		static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian) {
			jacobian.Redim(nodes.Rows(),0);
		}
		static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
		
	public:
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<long>& nodeindexes,
										  int matid, long& index);
	};
	
    template<class T>
    inline void TPZGeoPoint::X(const TPZFMatrix<REAL> &coord,TPZVec<T> &loc,TPZVec<T> &result){
        int i;
        for (i=0;i<coord.Rows();i++){
            result[i] = coord.GetVal(i,0);
        }
    }
    
    template<class T>
    inline void TPZGeoPoint::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        DebugStop();
    }
    

};

#endif

