/**
 * @file TPZGeoTriangle
 * @brief Contains the TPZGeoTriangle class which implements the geometry of a triangle element.
 */

#ifndef TPZGEOTRIANGLEH
#define TPZGEOTRIANGLEH

#include "pzvec.h"
#include "pzeltype.h"
#include "pznoderep.h"
#include "tpztriangle.h"
#include "pzfmatrix.h"

#include <string>
#include <map>

class TPZGeoEl;
class TPZGeoMesh;


namespace pzgeom {
	
	/**
	 * @ingroup geometry
	 * @brief Implements the geometry of a triangle element. \ref geometry "Geometry"
	 */
	class TPZGeoTriangle : public TPZNodeRep<3, pztopology::TPZTriangle> 
	{
	public:
        typedef pztopology::TPZTriangle Top;
		/** @brief Number of corner nodes */
		enum {NNodes = 3};
		
		/** @brief Constructor with list of nodes */
		TPZGeoTriangle(TPZVec<int64_t> &nodeindexes) : TPZRegisterClassId(&TPZGeoTriangle::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTriangle>(nodeindexes)
		{
		}
		
		/** @brief Empty constructor */
		TPZGeoTriangle() : TPZRegisterClassId(&TPZGeoTriangle::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTriangle>()
		{
		}
		
		/** @brief Constructor with node map */
		TPZGeoTriangle(const TPZGeoTriangle &cp,
					   std::map<int64_t,int64_t> & gl2lcNdMap) : TPZRegisterClassId(&TPZGeoTriangle::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp,gl2lcNdMap)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTriangle(const TPZGeoTriangle &cp) : TPZRegisterClassId(&TPZGeoTriangle::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}
		
		/** @brief Copy constructor */
		TPZGeoTriangle(const TPZGeoTriangle &cp, TPZGeoMesh &) : TPZRegisterClassId(&TPZGeoTriangle::ClassId),
        TPZNodeRep<NNodes,pztopology::TPZTriangle>(cp)
		{
		}


        void Jacobian(const TPZFMatrix<REAL> & coord, TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);

        static bool IsLinearMapping(int side)
        {
            return true;
        }

		/** @brief Returns the type name of the element */
		static std::string TypeName() { return "Triangle";}


        template<class T>
        static void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &x) {

            TPZFNMatrix<3,T> phi(3,1);
            TPZFNMatrix<6,T> dphi(2,3);
            TShape(loc,phi,dphi);
            int space = nodes.Rows();

            for(int i = 0; i < space; i++) {
                x[i] = 0.0;
                for(int j = 0; j < 3; j++) {
                    x[i] += phi(j,0)*nodes.GetVal(i,j);
                }
            }
        }


        /** @brief Compute gradient of x mapping from element nodes and local parametric coordinates */
        template<class T>
        static void GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);


//        /** @brief Compute the jacoabina associated to the x mapping from local parametric coordinates  */
//        static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
//                             TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);

		// /**
		//  * @brief Method which creates a geometric boundary condition
		//  * element based on the current geometric element, \n
		//  * a side and a boundary condition number
		//  */
		// static  TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig,int side,int bc);

        /** @brief Implementation of Hdiv space*/
        static	void ComputeNormal(TPZVec<REAL> &p1, TPZVec<REAL> &p2,TPZVec<REAL> &p3,TPZVec<REAL> &result);

        static	void VectorialProduct(TPZVec<REAL> &v1, TPZVec<REAL> &v2,TPZVec<REAL> &result);

//        static void VecHdiv(TPZFMatrix<REAL> & coord, TPZFMatrix<REAL> & fDeformedDirections,TPZVec<int> &sidevector);

        /* brief compute the vectors for defining an HDiv approximation space */
//        void VecHdiv(const TPZGeoEl &gel,TPZFMatrix<REAL> &NormalVec,TPZVec<int> & VectorSide) const
//        {
//            TPZFNMatrix<3*NNodes> coord(3,NNodes);
//            CornerCoordinates(gel, coord);
//            VecHdiv(coord,NormalVec,VectorSide);
//        }

        public:
            int ClassId() const override;
            void Read(TPZStream &buf, void *context) override;
            void Write(TPZStream &buf, int withclassid) const override;


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
		// 								  TPZVec<int64_t>& nodeindexes, int matid, int64_t& index);

	};


    template<class T>
    inline void TPZGeoTriangle::GradX(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){

        int space = nodes.Rows();
        int ncol = nodes.Cols();

        gradx.Resize(space,2);
        gradx.Zero();

#ifdef PZDEBUG
        if(/* nrow != 3 || */ ncol  != 3){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x3." << std::endl;
            DebugStop();
        }

#endif
        TPZFNMatrix<3,T> phi(NNodes,1);
        TPZFNMatrix<6,T> dphi(2,NNodes);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 3; i++)
        {
            for(int j = 0; j < space; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
                gradx(j,1) += nodes.GetVal(j,i)*dphi(1,i);

            }
        }

    }

};

#endif
