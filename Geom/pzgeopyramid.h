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
        typedef pztopology::TPZPyramid Top;
		/** @brief Number of corner nodes */
		enum {NNodes = 5};
		
            int ClassId() const override;
            void Read(TPZStream &buf, void *context) override;
            void Write(TPZStream &buf, int withclassid) const override;

                
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



        /** @brief Computes the jacobian*/
        static void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL>& par, TPZFMatrix<REAL> &jacobian, TPZFMatrix<REAL> &axes, REAL &detjac, TPZFMatrix<REAL> &jacinv);

        /** @brief Computes the geometric location*/
        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x);
        
        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);


		
		
	public:
        
        /// create an example element based on the topology
        /* @param gmesh mesh in which the element should be inserted
         @param matid material id of the element
         @param lowercorner (in/out) on input lower corner o the cube where the element should be created, on exit position of the next cube
         @param size (in) size of space where the element should be created
         */
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);

		/** @brief Creates a geometric element according to the type of the father element */
		// static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
		// 								  TPZVec<int64_t>& nodeindexes,
		// 								  int matid,
		// 								  int64_t& index);
	};

    
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
