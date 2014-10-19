//
//  TestHDivMesh.h
//  PZ
//
//  Created by Philippe Devloo on 10/19/14.
//
//

#ifndef __PZ__TestHDivMesh__
#define __PZ__TestHDivMesh__

#include <stdio.h>
#include "pzcompel.h"

int CompareShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

int CompareSideShapeFunctions(TPZCompElSide celsideA, TPZCompElSide celsideB);

void TestMesh(TPZCompMesh *cmesh);

#endif /* defined(__PZ__TestHDivMesh__) */
