#ifndef TPZARC3D_H
#define TPZARC3D_H

// #include "pzfmatrix.h"
// #include "pzvec.h"
// #include "pzgmesh.h"


//#include "pzgeoel.h"
#include "pznoderep.h"
// #include "pzgnode.h"
#include "tpzline.h"

#include <iostream>

using namespace std;
using namespace pzgeom;
using namespace pztopology;

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
		this->fICnBase = cp.fICnBase;
		this->fIBaseCn = cp.fIBaseCn;
		this->fCenter3D = cp.fCenter3D;
		this->fRadius = cp.fRadius;		
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
		int nnod = nodeindexes.NElements();
		if(nnod != 3)
		{
			std::cout << "Arc geometry created with " << nnod << " nodes, bailing out\n";
			DebugStop();
		}
		TPZFMatrix coord(3,nnod);
		int nod, co;
		for(nod=0; nod<3; nod++)
		{
			for(co=0; co<3; co++)
			{
				coord(co,nod) = mesh.NodeVec()[nodeindexes[nod]].Coord(co);
			}
		}
		ComputeAtributes(coord);
    }

    void X(TPZFMatrix &coord,TPZVec<REAL> &loc,TPZVec<REAL> &result);
    void Jacobian(TPZFMatrix &coord, TPZVec<REAL> &par, TPZFMatrix &jacobian, TPZFMatrix &axes, REAL &detjac, TPZFMatrix &jacinv);

    static std::string TypeName() { return "Linear";}
    static TPZGeoEl * CreateBCGeoEl(TPZGeoEl *orig, int side,int bc);
	
public:
	/**
	 * Creates a geometric element according to the type of the father element
	 */
	static TPZGeoEl *CreateGeoElement(TPZGeoMesh &mesh, MElementType type,
									  TPZVec<int>& nodeindexes,
									  int matid,
									  int& index);
	
	
protected:

    void ComputeAtributes(TPZFMatrix &coord);
    void ComputeR2Points(TPZFMatrix &coord, double &xa, double &ya, double &xb, double &yb, double &angle) const;
    double ArcAngle(TPZFMatrix &coord, double xa, double ya, double xb, double yb) const;


    /** Atributes */
    TPZFMatrix fICnBase;
    TPZFMatrix fIBaseCn;
    TPZVec< REAL > fCenter3D;
    double fRadius;
};



//#include "pzgeoelrefless.h.h"
//#include "tpzgeoelrefpattern.h.h"
//#include "pznoderep.h.h"

///CreateGeoElement -> TPZArc3D

#define TPZGEOELEMENTARC3DID 350
template<>
inline int TPZGeoElRefPattern<TPZArc3D>::ClassId() const {
	return TPZGEOELEMENTARC3DID;
}


//         template class pzgeom::TPZNodeRep<3,TPZArc3D>;
//template class TPZGeoElRefLess<TPZArc3D>;


#endif
