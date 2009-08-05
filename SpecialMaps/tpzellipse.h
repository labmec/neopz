#ifndef TPZELLIPSE_H
#define TPZELLIPSE_H

#include "pzgeoel.h"
#include "pznoderep.h"
#include "tpzline.h"

#include <iostream>

class TPZGeoMesh;

using namespace std;
using namespace pzgeom;
using namespace pztopology;

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

/** ESTA CLASSE SERÁ DELETADA!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! */
/** Será substituída pela classe TPZEllipse3D */

class TPZEllipse : public TPZNodeRep<2,TPZLine> {

private:

  void GetNodeCoord(TPZGeoMesh &mesh, TPZFMatrix &nodes);

public:

    enum {NNodes = 2};

    bool IsLinearMapping() const { return false; }

    TPZEllipse(const TPZEllipse &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
    }

    TPZEllipse() : TPZNodeRep<NNodes,pztopology::TPZLine>(){
    }

    TPZEllipse(const TPZEllipse &cp) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
    }

    TPZEllipse(const TPZEllipse &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
    }

    TPZEllipse(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes){
      this->AdjustNodeCoordinates(mesh);
    }

  /** Adjust node coordinate in case the non-linear mapping
   * changes node coordinates
   * It happens for instance in TPZEllipse
   */
  virtual void AdjustNodeCoordinates(TPZGeoMesh &mesh);

    void SetAxes(double a, double b);
    void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
    void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);

    static std::string TypeName() { return "Linear";}
    static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);

    /**
     * fA @param - Half of Ellipse Axe in X direction
     * fB @param - Half of Ellipse Axe in Y direction
     */
    double fA;
    double fB;
};

#endif
