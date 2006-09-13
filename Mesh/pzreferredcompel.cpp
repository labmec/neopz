// -*- c++ -*-

// $Id: pzreferredcompel.cpp,v 1.5 2006-09-13 19:18:25 cesar Exp $


#include "pzreferredcompel.h"
#include "pzelctemp.h"
#include "pzintel.h"
#include "TPZCompElDisc.h"

#include "pzquad.h"
#include "TPZGeoElement.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h" 
#include "TPZGeoLinear.h"
#include "pzshapelinear.h"
#include "pzgeoquad.h"
#include "pzshapequad.h"
#include "pzgeotriangle.h"
#include "pzshapetriang.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "pzgeoprism.h"
#include "pzshapeprism.h"
#include "pzgeopyramid.h"
#include "pzshapepiram.h"
#include "pzgeotetrahedra.h"
#include "pzshapetetra.h"

#include "pzmaterial.h"
#include "pzelmat.h"
#include "pzgeoel.h"
#include "pzcmesh.h"

#include "pzlog.h"



#include "pztempmat.h"

#ifdef DEBUG
  #define DEBUG2
#endif


#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.mesh.tpzcompel"));
#endif

using namespace std;

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::TPZReferredCompEl(TPZCompMesh &mesh, TPZGeoEl *gel, int &index):TCOMPEL(mesh, gel,index){

}//method

template<class TCOMPEL>
TPZReferredCompEl<TCOMPEL>::~TPZReferredCompEl(){

}//method

template< class TCOMPEL >
void TPZReferredCompEl< TCOMPEL >::ComputeSolution(TPZVec<REAL> &qsi, TPZFMatrix &phi, TPZFMatrix &dphix, TPZVec<REAL> &u, TPZFMatrix &du){

#ifdef DEBUG2
  TPZGeoEl * ref = this->Reference();
  if (!ref){
    std::stringstream mess; mess << __PRETTY_FUNCTION__ << " ERROR";
    LOGPZ_ERROR(logger, mess.str());
    PZError << mess.str() << std::endl;
  }
  TPZCompEl * othercel = ref->Reference();
//  MElementType tipo = othercel->Type();
  if (!othercel){
    std::stringstream mess; mess << __PRETTY_FUNCTION__ << " ERROR";
    LOGPZ_ERROR(logger, mess.str());
    PZError << mess.str() << std::endl;
  }  
#endif

  TPZManVector<REAL> ThisSol;
  TPZFNMatrix<100> ThisDSol(10,10);
  TCOMPEL::ComputeSolution(qsi, phi, dphix, ThisSol, ThisDSol);
  
  TPZManVector<REAL> OtherSol;
  TPZFNMatrix<100> OtherDSol(10,10);
  
  TPZCompEl * other = this->Reference()->Reference();
  TPZFMatrix axes (3,3,0.);

  other->ComputeSolution(qsi, OtherSol, OtherDSol,axes);

//  TPZManVector<REAL> AllSol(u.NElements() + OtherSol.NElements());
  u.Resize( ThisSol.NElements() + OtherSol.NElements() );
  for(int ii = 0; ii < ThisSol.NElements(); ii++){
    u[ii] = ThisSol[ii];
  }//for ii
  int ThisNEl = ThisSol.NElements();  
  for(int ii = 0; ii < OtherSol.NElements(); ii++){
    u[ii+ThisNEl] = OtherSol[ii];
  }//for ii
  
  const int dim = this->Reference()->Dimension();
  //TPZFNMatrix<100> AllDSol(dim, du.Cols() + OtherDSol.Cols());
  du.Redim(dim, ThisDSol.Cols() + OtherDSol.Cols());
  for(int ii = 0; ii < dim; ii++){
    for(int jj = 0; jj < ThisDSol.Cols(); jj++){
      du(ii,jj) = ThisDSol(ii,jj);
    }//for jj
  }//for ii
  int ThisCols = ThisDSol.Cols();
  for(int ii = 0; ii < dim; ii++){
    for(int jj = 0; jj < OtherDSol.Cols(); jj++){
      du(ii,ThisCols+jj) = OtherDSol(ii,jj);
    }//for jj
  }//for ii

}//method

using namespace pzshape;
using namespace pzgeom;

template class TPZReferredCompEl< TPZCompElDisc >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >;
template class TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >;

TPZCompEl * CreateReferredPointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPoint,TPZShapePoint> >(mesh,gel,index);
}

TPZCompEl * CreateReferredLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoLinear,TPZShapeLinear> >(mesh,gel,index);
}

TPZCompEl * CreateReferredQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoQuad,TPZShapeQuad> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTriangle,TPZShapeTriang> >(mesh,gel,index);
}

TPZCompEl * CreateReferredCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoCube,TPZShapeCube> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPrism,TPZShapePrism> >(mesh,gel,index);
}

TPZCompEl * CreateReferredPyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoPyramid,TPZShapePiram> >(mesh,gel,index);
}

TPZCompEl * CreateReferredTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZIntelGen<TPZGeoTetrahedra,TPZShapeTetra> >(mesh,gel,index);
}

TPZCompEl * CreateReferredDisc(TPZGeoEl *gel,TPZCompMesh &mesh,int &index) {
  return new TPZReferredCompEl< TPZCompElDisc >(mesh,gel,index);
}
