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
#include "pzfmatrix.h"

#include <string>

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
		/** @brief Number of corner nodes */
		enum {NNodes = 4};

                virtual int ClassId() const;
                
                void Read(TPZStream& buf, void* context);
                
                void Write(TPZStream& buf, int withclassid) const;

                
		/** @brief Constructor with list of nodes */
		TPZGeoTetrahedra(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoTetrahedra::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoTetrahedra() : TPZRegisterClassId(&TPZGeoTetrahedra::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTetrahedron>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp,
						 std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoTetrahedra::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp) : TPZRegisterClassId(&TPZGeoTetrahedra::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTetrahedra(const TPZGeoTetrahedra &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoTetrahedra::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTetrahedron>(cp)
		{
		}
		
        static bool IsLinearMapping(int side)
        {
            return true;
        }
        
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Tetrahedron";}
		
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
            TShape(loc, phi, dphi);
        }
		
        /* @brief Compute x mapping from local parametric coordinates */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &x) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,x);
        }
        
        /** @brief Compute gradient of x mapping from local parametric coordinates */
        template<class T>
        void GradX(const TPZGeoEl &gel, TPZVec<T> &loc, TPZFMatrix<T> &gradx) const
        {
            TPZFNMatrix<3*NNodes> coord(3,NNodes);
            CornerCoordinates(gel, coord);
            int nrow = coord.Rows();
            int ncol = coord.Cols();
            TPZFMatrix<T> nodes(nrow,ncol);
            for(int i = 0; i < nrow; i++)
            {
                for(int j = 0; j < ncol; j++)
                {
                    nodes(i,j) = coord(i,j);
                }
            }
            
            GradX(nodes,loc,gradx);
        }
        
        /** @brief Compute x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
        /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
        template<class T>
        static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
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
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

		/** @brief Creates a geometric element according to the type of the father element */
		static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
										  TPZVec<int64_t>& nodeindexes,
										  int matid,
										  int64_t& index);
	};
    
    template<class T>
    inline void TPZGeoTetrahedra::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T qsi = loc[0], eta = loc[1] , zeta  = loc[2];
        
        phi(0,0)  = 1.0-qsi-eta-zeta;
        phi(1,0)  = qsi;
        phi(2,0)  = eta;
        phi(3,0)  = zeta;
        
        dphi(0,0) = -1.0;
        dphi(1,0) = -1.0;
        dphi(2,0) = -1.0;
        dphi(0,1) =  1.0;
        dphi(1,1) =  0.0;
        dphi(2,1) =  0.0;
        dphi(0,2) =  0.0;
        dphi(1,2) =  1.0;
        dphi(2,2) =  0.0;
        dphi(0,3) =  0.0;
        dphi(1,3) =  0.0;
        dphi(2,3) =  1.0;
        
    }
    
    template<class T>
    inline void TPZGeoTetrahedra::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        TPZFNMatrix<4,T> phi(NNodes,1);
        TPZFNMatrix<12,T> dphi(3,NNodes);
        TShape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < NNodes; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
    }
    
    
    template<class T>
    inline void TPZGeoTetrahedra::GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 || ncol  != 4){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x4." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<4,T> phi(NNodes,1);
        TPZFNMatrix<12,T> dphi(3,NNodes);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 4; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
                gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);
                gradx(j,2) += nodes.GetVal(j,i)*dphi(2,i);
            }
        }
        
    }
	
};

#endif 
