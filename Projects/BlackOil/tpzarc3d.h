#ifndef TPZARC3D_H
#define TPZARC3D_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
// #include "pzgmesh.h"
#include "pzgeoel.h"
#include "pznoderep.h"
// #include "pzgnode.h"
#include "tpzline.h"

#include <iostream>

  /**
  / Class made by Paulo Cesar de Alvarenga Lucci (Caju)
  / LabMeC - FEC - UNICAMP
  / 2007
 */

class TPZArc3D : public TPZNodeRep<3,TPZLine> {

public:

    enum {NNodes = 3};

    bool IsLinearMapping() const { return false; }

    TPZArc3D(const TPZArc3D &cp,std::map<int,int> & gl2lcNdMap) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp,gl2lcNdMap){
    }

    TPZArc3D() : TPZNodeRep<NNodes,pztopology::TPZLine>(){
    }

    TPZArc3D(const TPZArc3D &cp) : TPZNodeRep<NNodes,pztopology::TPZLine>(cp){
          this->fICnBase = cp.fICnBase;
          this->fIBaseCn = cp.fIBaseCn;
          this->fCenter3D = cp.fCenter3D;
          this->fRadius = cp.fRadius;
    }

    TPZArc3D(const TPZArc3D &cp, TPZGeoMesh &) : TPZNodeRep<NNodes, pztopology::TPZLine>(cp){
          this->fICnBase  = cp.fICnBase;
          this->fIBaseCn  = cp.fIBaseCn;
          this->fCenter3D = cp.fCenter3D;
          this->fRadius   = cp.fRadius;
    }

    TPZArc3D(TPZVec<int> &nodeindexes, TPZGeoMesh &mesh) : TPZNodeRep<NNodes,pztopology::TPZLine>(nodeindexes){
    }

    void X(TPZFMatrix &nodes,TPZVec<REAL> &loc,TPZVec<REAL> &result);
    void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);

    static std::string TypeName() { return "Linear";}
    static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);

protected:

    void ComputeAtributes(TPZFMatrix &coord);
    void ComputeR2Points(TPZFMatrix &coord, double &xa, double &ya, double &xb, double &yb, double &angle);
    double ArcAngle(TPZFMatrix &coord, double xa, double ya, double xb, double yb);


    /** Atributes */
    TPZFMatrix fICnBase;
    TPZFMatrix fIBaseCn;
    TPZVec< REAL > fCenter3D;
    double fRadius;
};

#endif
