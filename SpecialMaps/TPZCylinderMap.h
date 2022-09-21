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
    
    template<class TGeo>
    class TPZCylinderMap : public TGeo
    {
        /// corner coordinates in cylindrical system (theta, z)
        TPZFNMatrix<TGeo::NNodes*3,REAL> fCornerCo;
        
        /// node around which we rotate the coordinates
        TPZManVector<REAL,3> fOrigin;
        
        /// radius of the cylinder
        REAL fRadius;
        
        /// rotation matrix
        TPZFNMatrix<9,REAL> fRotation;
        
    public:
        
        TPZCylinderMap() : TGeo(), fCornerCo(2*TGeo::NNodes,-1), fOrigin(3,0.), fRadius(0.), fRotation(3,3,0.)
        {
            fRotation.Identity();
        }
        
        TPZCylinderMap(TPZVec<int64_t> &nodeindices) : TGeo(nodeindices), fCornerCo(2*TGeo::NNodes,-1), fOrigin(3,0.), fRadius(0.), fRotation(3,3,0.)
        {
            fRotation.Identity();
        }
        
        TPZCylinderMap(const TPZCylinderMap &cp) : TGeo(cp), fCornerCo(cp.fCor), fRadius(cp.fRadius), fRotation(cp.fRotation)
        {
        }
        
        TPZCylinderMap &operator=(const TPZCylinderMap &cp)
        {
            TGeo::operator=(cp);
            fCornerCo = cp.fCornerCo;
            fOrigin = cp.fOrigin;
            fRadius = cp.fRadius;
            fRotation = cp.fRotation;
            return *this;
        }
        
        void SetOrigin(TPZVec<REAL> &origin, REAL radius)
        {
            fOrigin = origin;
            fRadius = radius;
        }

        /// axis direction with the vertical axis
        void SetCylinderAxis(const TPZVec<REAL> &axis);
        
        /// compute the corner coordinates of the corner nodes
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
        void GradX(const TPZGeoEl &gel, TPZVec<T> &par, TPZFMatrix<T> &gradx) const
        {
            
            
            TPZFNMatrix<6,T> DxDphi(3,2,0.), gradphi(2,2), gradxloc(3,2);
            TPZManVector<T,3> ft(2,0.);
            TGeo::X(fCornerCo,par,ft);
            TGeo::GradX(fCornerCo, par, gradphi);
            
            gradxloc(0,0) = -fRadius*sin(ft[0])*gradphi(0,0);
            gradxloc(1,0) = fRadius*cos(ft[0])*gradphi(0,0);
            gradxloc(2,0) = gradphi(1,0);
            gradxloc(0,1) = -fRadius*sin(ft[0])*gradphi(0,1);
            gradxloc(1,1) = fRadius*cos(ft[0])*gradphi(0,1);
            gradxloc(2,1) = gradphi(1,1);

            for (int64_t i=0; i<3; i++) {
                for (int64_t j=0; j<2; j++) {
                    gradx(i,j) = 0.;
                    for (int64_t k=0; k<3; k++) {
                        gradx(i,j) += fRotation.GetVal(k,i)*gradxloc(k,j);
                    }
                }
            }
        }
        
        template<class T>
        void X(const TPZFMatrix<REAL> &nodes,TPZVec<T> &loc,TPZVec<T> &result) const
        {
            TPZManVector<T,2> localcylinder(2);
            TGeo::X(this->fCornerCo,loc,localcylinder);
            TPZManVector<T,3> localcartesian(3);
            
            localcartesian[0] = fRadius*cos(localcylinder[0]);
            localcartesian[1] = fRadius*sin(localcylinder[0]);
            localcartesian[2] = localcylinder[1];
            
            for (int64_t i=0; i<3; i++) {
                result[i] = fOrigin[i];
                for (int64_t j=0; j<3; j++) {
                    result[i] += fRotation.GetVal(j,i)*localcartesian[j];
                }
            }
        }
        
        
        // static TPZGeoEl *CreateBCGeoEl(TPZGeoEl *gel, int side,int bc);
        
        /** @brief Creates a geometric element according to the type of the father element */
        // static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
        //                                   TPZVec<int64_t>& nodeindexes,
        //                                   int matid,
        //                                   int64_t& index);
        
        void Read(TPZStream& buf, void* context) override
        {
            TGeo::Read(buf,0);
            fCornerCo.Read(buf,0);
        }
        
        void Write(TPZStream &buf, int withclassid) const override
        {
            TGeo::Write(buf, withclassid);
            fCornerCo.Write(buf,0);
        }
        
        static void InsertExampleElement(TPZGeoMesh &gmesh, int matid, TPZVec<REAL> &lowercorner, TPZVec<REAL> &size);
        

    };
    
};

#endif /* TPZCylinder_hpp */
