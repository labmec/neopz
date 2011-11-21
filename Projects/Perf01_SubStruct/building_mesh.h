#ifndef UTIL_H

#include "pzvtkmesh.h"
#include <string.h>

void InsertElasticity(TPZAutoPointer<TPZCompMesh> mesh);

TPZGeoMesh* BuildBuildingMesh(const char* input_filename);

#endif
 
