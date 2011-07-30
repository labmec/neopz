/**
 * \file
 * @brief Creating TPZStack classe from template to geometric and computational elements.
 */
// $Id: pzstackmesh.cpp,v 1.1.1.1 2003-02-04 16:45:27 cantao Exp $

#include "pzstack.h"

#ifdef TPZEQNARRAY
/** Defining TPZEqnArray templates 15/08/2000 -> Longhin */
#include <tpzeqnarray.h>

template class TPZStack<TPZEqnArray>;
template class TPZStack<TPZEqnArray *>;

#endif // TPZEQNARRAY

#ifdef PZENVIRONMENT

#include "pzcompel.h"
#include "pzgeoel.h"
#include "pzelmat.h"

class TPZFMatrix;
template class TPZStack<TPZFMatrix *>;
class TPZMatPlaca2;
template class TPZStack<TPZMatPlaca2 *>;
class TPZCompEl;
template class TPZStack<TPZCompEl *>;
class TPZGeoEl;
template class TPZStack<TPZGeoEl *>;
class TPZCompCloneMesh;
template class TPZStack<TPZCompCloneMesh *>;

template class TPZStack<TPZCompElSide>;
template class TPZStack<TPZGeoElSide>;
class TPZCompMesh;
template class TPZStack<TPZCompMesh *>;
class TPZMatrixSolver;
template class TPZStack<TPZMatrixSolver *>;

//Including ElementMatrix
struct TPZElementMatrix;
template class TPZStack<TPZElementMatrix *>;

#endif // PZENVIRONMENT

#ifdef ELLIPS
#include "TElMeasurement.h"
template class TPZStack<TElMeasurement *>;
#endif // ELLIPS

/*
 #ifdef ARXDBX
 #include "dbid.h"
 template class TPZStack<AcDbObjectId>;
 #endif
 */

#ifdef ELLIPSARX
#include "TElRxPointerConv.h"
template class TPZStack<TElRxPointerConv::TG>;
#endif // ELLIPSARX

//--| PZ |----------------------------------------------------------------------
