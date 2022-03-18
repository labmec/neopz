//
//  TPZCylinder.cpp
//  pz
//
//  Created by Philippe Devloo on 04/05/18.
//

#include "TPZCylinder.h"
#include "pzvec_extras.h"
#include "tpzgeoelmapped.h"

namespace pzgeom {
    /// axis direction with the vertical axis
    template<class TGeo>
    void TPZCylinderMap<TGeo>::SetCylinderAxis(const TPZVec<REAL> &axis)
    {
        // build two normal vectors
        std::map<REAL,int> vals;
        for(int i=0; i<3; i++) vals[fabs(axis[i])] = i;
        TPZManVector<REAL,3> gen(3,0.), orto1(3,0.), orto2(3,0.);
        auto it = vals.rbegin();
        gen[it->second] = 1.;
        Cross(axis, gen, orto1);
        REAL normax = sqrt(Norm(axis));
        REAL normort1 = sqrt(Norm(orto1));
        for (int i=0; i<3; i++) {
            axis[i] /= normax;
            orto1[i] /= normort1;
        }
        Cross(orto1, axis, orto2);
        for (int i=0; i<3; i++) {
            fRotation(0,i) = orto1[i];
            fRotation(1,i) = orto2[i];
            fRotation(2,i) = axis[i];
        }
    }
    
    /// compute the corner coordinates of the corner nodes
    template<class TGeo>
    void TPZCylinderMap<TGeo>::ComputeCornerCoordinates(TPZGeoMesh &gmesh)
    {
        int nnodes = TGeo::NNodes;
        for (int in=0; in<nnodes; in++) {
            int64_t nodeindex = TGeo::fNodeIndexes[in];
            TPZManVector<REAL,3> co(3);
            gmesh.NodeVec()[nodeindex].GetCoordinates(co);
            TPZManVector<REAL,3> localco(3);
            for (int i=0; i<3; i++) {
                localco[i] = 0.;
                for (int j=0; j<3; j++) {
                    localco[i] += fRotation(i,j)*co[i]-fOrigin[i];
                }
            }
            REAL radius = sqrt(localco[0]*localco[0]+localco[1]*localco[1]);
            if(fabs(radius-fRadius) > 1.e-6) DebugStop();
            REAL theta = atan2(localco[1],localco[0]);
            REAL z = localco[2];
            fCornerCo(0,in) = theta;
            fCornerCo(1,in) = z;
        }
    }
    
    
    

};
