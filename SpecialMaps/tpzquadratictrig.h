/**
 * @file
 * @brief Contains the TPZQuadraticTrig class which defines a triangular geometric element with quadratic map.
 */
 
#ifndef TPZQUADRATICTRIG_H
#define TPZQUADRATICTRIG_H

#include "pznoderep.h"
#include "tpztriangle.h"

class TPZGeoMesh;

namespace pzgeom {
    
	/**
	 * @author Paulo Cesar de Alvarenga Lucci (Caju)
	 * @since 2007
	 * @ingroup geometry
	 * @brief Defines a triangular geometric element with quadratic map. \ref geometry "Geometry"
	 */
    class TPZQuadraticTrig : public pzgeom::TPZNodeRep<6,pztopology::TPZTriangle> {
		
	public:
        typedef pztopology::TPZTriangle Top;
		/** @brief Number of nodes */
		enum {NNodes = 6};
                
                public:
int ClassId() const override;

        
        //irtual void ParametricDomainNodeCoord(int node, TPZVec<REAL> &nodeCoord);
        
		/** @brief It is a quadratic mapping */
        static bool IsLinearMapping(int side)
        {
            return false;
        }
        
		/** @brief Constructor for node indexes given */
		TPZQuadraticTrig(TPZVec<int64_t> &nodeindexes) : 
        TPZRegisterClassId(&TPZQuadraticTrig::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes) {
		}
		/** @brief Default constructor */
		TPZQuadraticTrig() : TPZRegisterClassId(&TPZQuadraticTrig::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>() {
		}
		/** @brief Copy constructor for node map given */
		TPZQuadraticTrig(const TPZQuadraticTrig &cp,std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZQuadraticTrig::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap) {
		}
		/** @brief Copy constructor */
		TPZQuadraticTrig(const TPZQuadraticTrig &cp) : TPZRegisterClassId(&TPZQuadraticTrig::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
		}
		/** @brief Copy constructor */
		TPZQuadraticTrig(const TPZQuadraticTrig &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZQuadraticTrig::ClassId),
        pzgeom::TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp) {
		}
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "QuadraticTriangle";}
		
		// /**
		//  * @brief Method which creates a geometric boundary condition 
		//  * element based on the current geometric element, 
		//  * a side and a boundary condition number
		//  */
		// static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);
		
        /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
        
		
		/* brief compute the coordinate of a point given in parameter space */
//        template<class T>
//        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            X(coord,loc,result);
//        }
        
        /** @brief Compute gradient of x mapping from local parametric coordinates */
//        template<class T>
//        void GradX(const TPZGeoEl &gel, TPZVec<T> &loc, TPZFMatrix<T> &gradx) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            GradX(coord,loc,gradx);
//        }
		
        template<class T>
        static void TShape(const TPZVec<T> &param,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
        template<class T>
		static void X(const TPZFMatrix<REAL> &coord, TPZVec<T> &par, TPZVec< T > &result);
        
        /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
		
	public:
		// /** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid, int64_t& index);

        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

	};

};

#endif
