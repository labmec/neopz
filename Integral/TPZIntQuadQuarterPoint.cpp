//
//  TPZIntQuadQuarterPoint.cpp
//  PZ
//
//  Created by Philippe Devloo on 6/1/15.
//
//

#include "TPZIntQuadQuarterPoint.h"


/// return the weight and position of the integration point
void TPZIntQuadQuarterPoint::Point(int ip, TPZVec<REAL> &pos, REAL &w) const
{
    TPZManVector<REAL,2> ptr(2);
    TPZIntQuad::Point(ip/2, ptr, w);
    REAL theta = (ptr[0]+1.)*M_PI/8.;
    REAL rrange = 1./cos(theta);
    REAL r = (ptr[1]+1.)*rrange/2.;
    w *= r*rrange*M_PI/4.;
    //std::cout << "ip " << ip << " ptr " << ptr << " theta " << theta << " rrange " << rrange << " r " << r << std::endl;
    
    REAL xi,eta;
    xi = r*cos(theta);
    eta = r*sin(theta);
    //std::cout << "xi " << xi << " eta " << eta << std::endl;
    if (ip%2) {
        pos[0] = -1.+2.*xi;
        pos[1] = -1.+2.*eta;
    }
    else
    {
        pos[0] = -1.+2.*eta;
        pos[1] = -1.+2.*xi;
    }
    switch (fCorner) {
        case 0:
            break;
        case 1:
            pos[0] = -pos[0];
            break;
        case 2:
            pos[0] = -pos[0];
            pos[1] = -pos[1];
            break;
        case 3:
            pos[1] = -pos[1];
            break;
        default:
            DebugStop();
            break;
    }
    //std::cout << "pos " << pos << std::endl;
}

/// set the integration order
void TPZIntQuadQuarterPoint::SetOrder(TPZVec<int> &ord,int type)
{
    TPZIntQuad::SetOrder(ord, type);
    fNPoints = TPZIntQuad::NPoints()*2;
}

