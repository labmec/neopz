#ifndef vigaH
#define vigaH

#include "TPZPullOutTestAnalysis.h"
//#include "TImportMesh.h"
//#include "TPZAxisymmetricMat.h"
#include "pzelast3dGD.h"
//#include "TPZMohrCoulombModel.h"
#include "pzvisualmatrix.h"
#include "pzskylstrmatrix.h"
#include "pzanalysis.h"
#include "pzstepsolver.h"
#include "TSWXGraphMesh.h"
#include "TSWXGraphElement.h"
#include "pzgeoel.h"
#include "tpzautopointer.h"
#include "TPZMaterial.h"
#include "pzbndcond.h"
#include "TPZRefPatternDataBase.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZFrontStructMatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzfstrmatrix.h"
#include "TPZCompElDisc.h"
#include "TPZInterfaceEl.h"
#include "pzintel.h"
#include "pznonlinanalysis.h"
#include "pzeltype.h"
//#include "Dialogs.hpp"


TPZGeoMesh * GeoMesh();

TPZCompMesh * CompMesh(TPZGeoMesh * gmesh, int pOrder);

void IncidHexa(int nx, int ny, int iel, TPZVec<int> &incid);
int CreateMidNode(int nodeA, int nodeB, TPZGeoMesh * gmesh, REAL Yoffset);




#endif
