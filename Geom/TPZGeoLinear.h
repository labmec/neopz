
#ifndef TPZGEOLINEARH
#define TPZGEOLINEARH

/**
 * @file TPZGeoLinear
 * @brief Implements the geometry of a one dimensional linear element.
 */

#include "pznoderep.h"

#include "pzvec.h"
#include "pzeltype.h"
#include "pzfmatrix.h"
#include "tpzline.h"

#include <string>


template<class TVar>
class TPZFMatrix;
class TPZGeoEl;
class TPZGeoMesh;

namespace pzgeom
{

	
    /**
     * @ingroup geometry
     * @brief Implements the geometry of a one dimensional linear element. \ref geometry "Geometry"
     */
    class TPZGeoLinear : public TPZNodeRep<2, pztopology::TPZLine> {
            
        public:
        
            /** @brief Number of corner nodes */
            enum {NNodes = 2};
            
            /** @brief Constructor with list of nodes */
            TPZGeoLinear(TPZVec<long> &nodeindexes) : TPZNodeRep<NNodes, pztopology::TPZLine>(nodeindexes)
            {
            }
            
            /** @brief Empty constructor */
            TPZGeoLinear() : TPZNodeRep<NNodes, pztopology::TPZLine>()
            {
            }
            
            /** @brief Constructor with node map */
            TPZGeoLinear(const TPZGeoLinear &cp, std::map<long,long> & gl2lcNdMap) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp,gl2lcNdMap)
            {
            }
            
            /** @brief Copy constructor */
            TPZGeoLinear(const TPZGeoLinear &cp) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
            {
            }
            
            /** @brief Copy constructor */
            TPZGeoLinear(const TPZGeoLinear &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp)
            {
            }
            
            /** @brief answer if the element side is a linear map */
            static bool IsLinearMapping(int side)
            {
                return true;
            }
            
            /** @brief Returns the type name of the element */
            static std::string TypeName() { return "Linear";}
        
            /** @brief Compute the shape being used to construct the X mapping from local parametric coordinates  */
            static void Shape(TPZVec<REAL> &loc,TPZFMatrix<REAL> &phi,TPZFMatrix<REAL> &dphi){
                TShape(loc, phi, dphi);
            }
        
            /* @brief Compute X mapping from local parametric coordinates */
            void X(const TPZGeoEl &gel,TPZVec<REAL> &loc,TPZVec<REAL> &x) const
            {
                TPZFNMatrix<3*NNodes> coord(3,NNodes);
                CornerCoordinates(gel, coord);
                X(coord,loc,x);
            }
        
            /** @brief Compute gradient of X mapping from local parametric coordinates */
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
            
            /* @brief Computes the jacobian of the map between the master element and deformed element */
            void Jacobian(const TPZGeoEl &gel,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv) const
            {
                TPZFNMatrix<3*NNodes> coord(3,NNodes);
                CornerCoordinates(gel, coord);
                Jacobian(coord, param, jacobian, axes, detjac, jacinv);
            }
        
            /** @brief Compute X mapping from element nodes and local parametric coordinates */
            static void X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &x);
        
            /** @brief Compute gradient of X mapping from element nodes and local parametric coordinates */
            template<class T>
            static void GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx);
        
            /** @brief Compute the shape being used to construct the x mapping from local parametric coordinates  */
            template<class T>
            static void TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi);
        
            /** @brief Compute the jacoabina associated to the x mapping from local parametric coordinates  */
            static void Jacobian(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &param,TPZFMatrix<REAL> &jacobian,
                                 TPZFMatrix<REAL> &axes,REAL &detjac,TPZFMatrix<REAL> &jacinv);
        
            /**
             * @brief Method which creates a geometric boundary condition
             * element based on the current geometric element, \n
             * a side and a boundary condition number
             */
            static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
            
            
        public:
        
        /** @brief Creates a geometric element according to the type of the father element */
        static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh,
                                          MElementType type,
                                          TPZVec<long>& nodeindexes,
                                          int matid,
                                          long& index);
        
    };
    
    template<class T>
    inline void TPZGeoLinear::TShape(TPZVec<T> &loc,TPZFMatrix<T> &phi,TPZFMatrix<T> &dphi) {
        T x = loc[0];
        phi(0,0) = (1.0-x)/2.;
        phi(1,0) = (1.0+x)/2.;
        dphi(0,0) = -0.5;
        dphi(0,1) = 0.5;
    }
    
    inline void TPZGeoLinear::X(const TPZFMatrix<REAL> &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &x){
        
        REAL xi = loc[0];
        int nrow = nodes.Rows();
        for(int i = 0; i < nrow; i++)
        {
            x[i] = nodes.GetVal(i,0)*(1.-xi)*0.5+nodes.GetVal(i,1)*(1.+xi)*0.5;
        }
        
    }
    
    
    template<class T>
    inline void TPZGeoLinear::GradX(const TPZFMatrix<T> &nodes,TPZVec<T> &loc, TPZFMatrix<T> &gradx){
        
        gradx.Resize(3,1);
        gradx.Zero();
        int nrow = nodes.Rows();
        int ncol = nodes.Cols();
#ifdef PZDEBUG
        if(nrow != 3 && ncol  != 2){
            std::cout << "Objects of incompatible lengths, gradient cannot be computed." << std::endl;
            std::cout << "nodes matrix must be 3x2." << std::endl;
            DebugStop();
        }
        
#endif
        TPZFNMatrix<3,T> phi(2,1);
        TPZFNMatrix<6,T> dphi(2,2);
        TShape(loc,phi,dphi);
        for(int i = 0; i < 2; i++)
        {
            for(int j = 0; j < 3; j++)
            {
                gradx(j,0) += nodes.GetVal(j,i)*dphi(0,i);
            }
        }
        
    }
	
};

#endif
