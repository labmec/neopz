
#include "postprocess.h"
#include "pzblackoil2p3d.h"
#include "pzblackoilanalysis.h"

#include "TPZGeoElement.h"
#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "TPZInterfaceEl.h"

#include "pzcompel.h"
#include "TPZCompElDisc.h"
#include "pzbndcond.h"


double Pressao(TPZBlackOilAnalysis &an, int matid){
  an.LoadSolution(an.Solution());
  TPZCompMesh * cmesh = an.Mesh();
  const int nel = cmesh->NElements();
  TPZVec<REAL> qsi(3), sol(1);
  double press = 0.;
  double AccVol = 0.;
  double locVol = 0.;
  for(int iel = 0; iel < nel; iel++){
    TPZCompEl * cel = cmesh->ElementVec()[iel];
    if(!cel) continue;
    if(cel->Material()->Id() != matid) continue;
    TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*>(cel);
    if(!face) continue;
    face->LeftElement()->Solution(qsi, TPZBlackOil2P3D::EOilPressure, sol);
    locVol = face->LeftElement()->Reference()->Volume();
    press += sol[0] * locVol;
    AccVol += locVol;
  }///iel
  double result = press/AccVol;
  return result;
}///method

void Vazao(TPZBlackOilAnalysis &an, int matid, double & VazaoAgua, double  & VazaoOleo){

}///method

