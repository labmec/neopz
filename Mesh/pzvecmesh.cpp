/**
 * @file
 * @brief Creates the TPZVec classes for computational and geometric structures.
 */
//$Id: pzvecmesh.cpp,v 1.2 2003-11-05 16:02:21 tiago Exp $

#include "pzvec.h"

#ifdef PZENVIRONMENT

class TPZCompEl;
class TPZConnect;
struct TPZConnectBC;
class TPZMatPlaca2;
class TPZFMatrix;
class TPZCompElQ2d;
template class TPZVec<TPZFMatrix<REAL> *>;
template class TPZVec<TPZCompEl *>;
template class TPZVec<TPZCompEl **>;
template class TPZVec<TPZCompElQ2d *>;
template class TPZVec<TPZConnect *>;
template class TPZVec<TPZConnectBC *>;
template class TPZVec<TPZMatPlaca2 *>;

class TPZGeoEl;
template class TPZVec<TPZGeoEl **>;
template class TPZVec<TPZGeoEl *>;
class TPZGeoNode;
template class TPZVec<TPZGeoNode *>;
template class TPZVec<TPZGeoNode **>;
struct TPZGeoNodeBC;
template class TPZVec<TPZGeoNodeBC *>;
class TPZMaterial;
template class TPZVec<TPZMaterial **>;
class TPZCosys;
template class TPZVec<TPZCosys *>;
template class TPZVec<TPZCosys **>;
struct TPZGeoElBC;
template class TPZVec<TPZGeoElBC *>;
class TPZBndCond;
template class TPZVec<TPZBndCond **>;

#include <pzgeoel.h>
template class TPZVec<TPZGeoElSide>;

#include "pzcompel.h"
template class TPZVec<TPZCompElSide>;

class TPZGraphEl;
template class TPZVec<TPZGraphEl ***>;
template class TPZVec<TPZGraphEl **>;
template class TPZVec<TPZGraphEl *>;
class TPZGraphNode;
//#include "pzgraphnode.h"
template class TPZVec<TPZGraphNode **>;
template class TPZVec<TPZGraphNode *>;
class TPZCompMesh;
template class TPZVec<TPZCompMesh *>;
class TPZMatrixSolver;
template class TPZVec<TPZMatrixSolver *>;

class TPZCompCloneMesh;
template class TPZVec<TPZCompCloneMesh *>;


#ifdef TPZEQNARRAY
/**Defining TPZEqnArray templates  -> 15/08/2000 by Longhin*/
#include <tpzeqnarray.h>
template class TPZVec<TPZEqnArray>;
template class TPZVec<TPZEqnArray *>;
#endif
#include "pzelmat.h"
struct TPZElementMatrix;
template class TPZVec<TPZElementMatrix *>;


#endif
