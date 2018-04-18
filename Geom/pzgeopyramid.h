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
#include "pzfmatrix.h"

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
		
            virtual int ClassId() const;
            void Read(TPZStream& buf, void* context);
            void Write(TPZStream& buf, int withclassid) const;

                
		/** @brief Constructor with list of nodes */
		TPZGeoPyramid(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoPyramid::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPyramid>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoPyramid() : TPZRegisterClassId(&TPZGeoPyramid::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPyramid>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoPyramid(const TPZGeoPyramid &cp,
					  std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoPyramid::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPyramid(const TPZGeoPyramid &cp) : TPZRegisterClassId(&TPZGeoPyramid::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoPyramid(const TPZGeoPyramid &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoPyramid::ClassId),
        TPZNodeRep<NNodes, pztopology::TPZPyramid>(cp)
		{
		}
        
        static bool IsLinearMapping(int side)
        {
            return true;
        }
		
		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Pyramid";}
		
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
//            int nrow = coord.Rows();
//            int ncol = coord.Cols();
//            TPZFMatrix<T> nodes(nrow,ncol);
//            for(int i = 0; i < nrow; i++)
//            {
//                for(int j = 0; j < ncol; j++)
//                {
//                    nodes(i,j) = coord(i,j);
//                }
//            }
            
            GradX(coord,loc,gradx);
        }

        /* @brief compute the jacobian of the map between the master element and deformed element */
        void Jacobian(const TPZGeoEl &gel, TPZVec<REAL> &param, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv) const {
            TPZFNMatrix<3 * NNodes> coord(3, NNodes);
            CornerCoordinates(gel, coord);
            Jacobian(coord, param, jacobian, axes, detjac, jacinv);
        }

        /** @brief Computes the jacobian*/
        static void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);

        /** @brief Computes the geometric location*/
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
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
    inline void TPZGeoPyramid::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T xi = loc[0], eta = loc[1] , zeta  = loc[2];
        
        if (zeta> 1.) {
            DebugStop();
        }
        
        T T0xz = .5*(1.-zeta-xi) / (1.-zeta);
        T T0yz = .5*(1.-zeta-eta) / (1.-zeta);
        T T1xz = .5*(1.-zeta+xi) / (1.-zeta);
        T T1yz = .5*(1.-zeta+eta) / (1.-zeta);
        if (IsZero(xi)) {
            T0xz = 0.5;
            T1xz = 0.5;
        }
        if (IsZero(eta)) {
            T0yz = 0.5;
            T1yz = 0.5;
        }
        T lmez = (1.-zeta);
        
        phi(0,0)  = T0xz*T0yz*lmez;
        phi(1,0)  = T1xz*T0yz*lmez;
        phi(2,0)  = T1xz*T1yz*lmez;
        phi(3,0)  = T0xz*T1yz*lmez;
        phi(4,0)  = zeta;
        
        T lmexmez = 1.-xi-zeta;
        T lmeymez = 1.-eta-zeta;
        T lmaxmez = 1.+xi-zeta;
        T lmaymez = 1.+eta-zeta;
        
        if (IsZero(lmez) && !IsZero(lmexmez) && !IsZero(lmeymez) &&
            !IsZero(lmaxmez) && !IsZero(lmaymez)) {
            DebugStop();
        }
        if (IsZero(lmez)) {
            lmexmez = 0.999;
            lmeymez = 0.999;
            lmaxmez = 0.999;
            lmaymez = 0.999;
            lmez = 0.001;
        }
        
        dphi(0,0) = -.25*lmeymez / lmez;
        dphi(1,0) = -.25*lmexmez / lmez;
        dphi(2,0) = -.25*(lmeymez+lmexmez-lmexmez*lmeymez/lmez) / lmez;
        dphi(0,1) =  .25*lmeymez / lmez;
        dphi(1,1) = -.25*lmaxmez / lmez;
        dphi(2,1) = -.25*(lmeymez+lmaxmez-lmaxmez*lmeymez/lmez) / lmez;
        dphi(0,2) =  .25*lmaymez / lmez;
        dphi(1,2) =  .25*lmaxmez / lmez;
        dphi(2,2) = -.25*(lmaymez+lmaxmez-lmaxmez*lmaymez/lmez) / lmez;
        dphi(0,3) = -.25*lmaymez / lmez;
        dphi(1,3) =  .25*lmexmez / lmez;
        dphi(2,3) = -.25*(lmaymez+lmexmez-lmexmez*lmaymez/lmez) / lmez;
        dphi(0,4) =  0.0;
        dphi(1,4) =  0.0;
        dphi(2,4) =  1.0;
        
    }
    
    template<class T>
    inline void TPZGeoPyramid::X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x){
        
        TPZFNMatrix<5,T> phi(5,1);
        TPZFNMatrix<15,T> dphi(3,5);
        TShape(loc,phi,dphi);
        int space = nodes.Rows();
        
        for(int i = 0; i < space; i++) {
            x[i] = 0.0;
            for(int j = 0; j < 5; j++) {
                x[i] += phi(j,0)*nodes.GetVal(i,j);
            }
        }
    }
    
    
    template<class T>
    inline void TPZGeoPyramid::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,3);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 || ncol  != 5){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x5." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<5,T> phi(5,1);
        TPZFNMatrix<15,T> dphi(3,5);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 5; i++)
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
