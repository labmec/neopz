/**
 * @file
 * @brief Contains the TPZReadMeshHR class which reads a mesh in a "human readable" format.
 */
/*****************************************************************************
 * O contedo desse arquivo �de propriedade do LabMeC-DES-FEC-UNICAMP e do
 * CENPES-Petrobras. 
 * O uso de qualquer parte ou do todo est�condicionado �expressa autoriza�o
 * dos propriet�ios.
 *****************************************************************************/

#ifndef PZREADMESHHR_H
#define PZREADMESHHR_H

#include "pzmanvector.h"
#include <string>

class TPZGeoMesh;
class TPZGeoEl;
class TPZCompMesh;
class TPZGeoElSide;

/**
 * @ingroup pre
 * @brief Reads a mesh in a "human readable" format, i.e. in text format and with coments. \ref pre "Getting Data"
 * @author Edimar Cesar Rylo
 * @since September 2006
 */
/** 
 * The lines that contains comments must start with a ":"
 * Note that this parser provides interface for read only 2D elasticity materials!
 */
class TPZHyperPlane
{
    
    /// center of the hyperplane
    TPZManVector<REAL,3> fCenter;
    
    /// normal to the hyperplane
    TPZManVector<REAL,3> fNormal;
    
public:

    TPZHyperPlane(TPZVec<REAL> &center, TPZVec<REAL> &normal) : fCenter(center), fNormal(normal)
    {
        
    }

    bool IsLeft(TPZVec<REAL> &point) const
    {
        REAL val=0.;
        for (int i=0; i<3; i++) {
            val += (point[i]-fCenter[i])*fNormal[i];
        }
        return val <= 0;
    }
    
    REAL Distance(const TPZVec<REAL> &point, TPZVec<REAL> &jacobian) const
    {
        REAL distance = 0;
        for (int i=0; i<3; i++) {
            distance += (point[i]-fCenter[i])*fNormal[i];
            jacobian[i] = fNormal[i];
        }
        return distance;
    }
};

class TPZHyperPlaneIntersect
{
public:
    /// Compute the intersection of the geometric mesh with the plane (only leaf elements)
    void Intersect(TPZGeoMesh &input, const TPZHyperPlane &plane, TPZGeoMesh &target);
    
    /// Compute the intersection between the edge and the plane
    REAL EdgeIntersect(const TPZGeoElSide &gelside, const TPZHyperPlane &plane);
    
private:
    /// reorder the nodes to form a convex figure
    void Reorder(TPZGeoEl *gel, TPZGeoMesh &target, TPZVec<std::pair<int,int> > &sidenodepair);
    /// a version which allows for more than 4 nodes
    int ReorderGeneral(TPZGeoMesh &target, TPZVec<std::pair<int,int> > &sidenodepair);
};

#endif
