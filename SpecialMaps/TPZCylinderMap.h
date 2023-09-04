//
//  TPZCylinder.hpp
//  pz
//
//  Created by Philippe Devloo on 04/05/18.
//

#ifndef TPZCylinder_hpp
#define TPZCylinder_hpp

#include <stdio.h>

#include "pzgeotriangle.h"

namespace pzgeom {
    /**
     @brief Implements elements in a cylindrical shell.
    It is based on storing information in cylindrical coordinates
    and adjusting to the mapped cylinder via translation + rotation.
    The rotation is computed based on the transformation between
    the reference cylinder axis ({0,0,1}) and the mapped cylinder.*/
    template<class TGeo>
    class TPZCylinderMap : public TGeo
    {
        /// corner coordinates in cylindrical system (theta, z)
        TPZFNMatrix<TGeo::NNodes*3,REAL> fCylindricalCo;
        
        /// node around which we rotate the coordinates
        TPZManVector<REAL,3> fOrigin;
        
        
        /**
           @brief Computes the rotation matrix of the cylinder.
           The reference cylinder has axis {0,0,1}. 
           This matrix is such that 
              fRotation . {0,0,1}
           equals to the  axis of the deformed cylinder.
         */
        TPZFNMatrix<9,REAL> fRotation;
        
    public:
        
        TPZCylinderMap() : TGeo(), fCylindricalCo(3,TGeo::NNodes,-1), fOrigin(3,0.), fRotation(3,3,0.)
        {
            fRotation.Identity();
        }
        
        TPZCylinderMap(TPZVec<int64_t> &nodeindices) : TGeo(nodeindices), fCylindricalCo(3,TGeo::NNodes,-1), fOrigin(3,0.), fRotation(3,3,0.)
        {
            fRotation.Identity();
        }
        
        TPZCylinderMap(const TPZCylinderMap &cp) : TGeo(cp), fCylindricalCo(cp.fCylindricalCo), fRotation(cp.fRotation)
        {
        }

        /** @brief Copy constructor for another mesh*/
        TPZCylinderMap(const TPZCylinderMap &cp,
                       TPZGeoMesh &destmesh) : TGeo(cp,destmesh), fCylindricalCo(cp.fCylindricalCo),
                                               fOrigin(cp.fOrigin),
                                               fRotation(cp.fRotation)
        {}

        /** @brief Copy constructor with node map*/
        TPZCylinderMap(const TPZCylinderMap &cp,
                       std::map<int64_t,int64_t> & gl2lcNdMap) : TGeo(cp,gl2lcNdMap), fCylindricalCo(cp.fCylindricalCo),
                                                                 fOrigin(cp.fOrigin),
                                                                 fRotation(cp.fRotation)
        {}
        
        TPZCylinderMap &operator=(const TPZCylinderMap &cp)
        {
            TGeo::operator=(cp);
            fCylindricalCo = cp.fCylindricalCo;
            fOrigin = cp.fOrigin;
            fRotation = cp.fRotation;
            return *this;
        }

        static bool IsLinearMapping(int side)
        {
            return false;
        }

        /// Sets center of the cylinder's base and radius
        void SetOrigin(const TPZVec<REAL> &origin)
        {
            fOrigin = origin;
        }

        /// axis direction with the vertical axis
        /// the last line indicates the z-direction
        void SetCylinderAxis(const TPZFMatrix<REAL> &axis);
        
        /// compute the cylindrical coordinates of the corner nodes
        void ComputeCornerCoordinates(TPZGeoMesh &gmesh);
        
        /** @brief Returns the type name of the element */
        static std::string TypeName() { return "CylinderMap";}
        
        /* @brief Computes the coordinate of a point given in parameter space */
        template<class T>
        void X(const TPZGeoEl &gel,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZFNMatrix<3*TGeo::NNodes> coord(3,TGeo::NNodes);
            CornerCoordinates(gel, coord);
            X(coord,loc,result);
        }
 
        template<class T>
        void ComputeGradRX(const TPZVec<T> &ftz, TPZFMatrix<T> &gradrx) const {
            gradrx(0,0) = cos(ftz[1]);
            gradrx(0,1) = -ftz[0]*sin(ftz[1]);
            gradrx(0,2) = ftz[0]*0.;
            gradrx(1,0) = sin(ftz[1]);
            gradrx(1,1) = ftz[0]*cos(ftz[1]);
            gradrx(1,2) = ftz[0]*0.;
            gradrx(2,0) = ftz[0]*0.;
            gradrx(2,1) = ftz[0]*0.;
            gradrx(2,2) = ftz[0]*0.+1.;
        }

        template<class T>
        void GradX(const TPZFMatrix<REAL> &nodes, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            
            
            TPZFNMatrix<9,T> gradphi(3,3), gradxloc(3,3);
            TPZManVector<T,3> rtz(3,0.);
            TGeo::X(fCylindricalCo,par,rtz);
            TGeo::GradX(fCylindricalCo, par, gradphi);
            ComputeGradRX(rtz,gradxloc);
            // melhorar
            /*
            gradxloc(0,0) = -fRadius*sin(ft[0])*gradphi(0,0);
            gradxloc(1,0) = fRadius*cos(ft[0])*gradphi(0,0);
            gradxloc(2,0) = gradphi(1,0);
            gradxloc(0,1) = -fRadius*sin(ft[0])*gradphi(0,1);
            gradxloc(1,1) = fRadius*cos(ft[0])*gradphi(0,1);
            gradxloc(2,1) = gradphi(1,1);
*/
            for (int64_t i=0; i<3; i++) {
                for (int64_t j=0; j<TGeo::Dimension; j++) {
                    gradx(i,j) = 0.;
                    for (int64_t k=0; k<3; k++) {
                        for (int l = 0; l<3; l++) {
                            gradx(i,j) += fRotation.GetVal(i,k)*gradxloc(k,l)*gradphi(l,j);
                        }
                    }
                }
            }
        }
        
        template<class T>
        void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZManVector<T,3> localcylinder(3);
            TGeo::X(this->fCylindricalCo,loc,localcylinder);
            TPZManVector<T,3> localcartesian(3);
            T Radius = localcylinder[0];
            T theta = localcylinder[1];
            localcartesian[0] = Radius*cos(theta);
            localcartesian[1] = Radius*sin(theta);
            localcartesian[2] = localcylinder[2];
            
            for (int64_t i=0; i<3; i++) {
                result[i] = fOrigin[i];
                for (int64_t j=0; j<3; j++) {
                    result[i] += fRotation.GetVal(i,j)*localcartesian[j];
                }
            }
        }
        
        
        // static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
        
        /** @brief Creates a geometric element according to the type of the father element */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t>& nodeindexes,
        //                                   int matid,
        //                                   int64_t& index);
        
        void Read(TPZStream& buf, void* context) override;
        
        void Write(TPZStream &buf, int withclassid) const override;
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

    };
    
};

#endif /* TPZCylinder_hpp */
