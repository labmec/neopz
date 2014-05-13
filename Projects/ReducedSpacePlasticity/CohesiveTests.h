//
//  TPZCohesiveTests.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef CohesiveTests
#define CohesiveTests

#include "pzreal.h"

class TPZGeoMesh;
class TPZCompMesh;
class TPZNonLinearAnalysis;

void ElastTest();
void ElastNLTest();
void ElastNLTestWithCohesive();
void CohesiveTwoLoads();
void OpenFractureTest();

void TestContinuousDiscontinuous();

TPZGeoMesh* CreateGeoMesh();
TPZGeoMesh* CreateGeoMeshCohe();
TPZGeoMesh* CreateGeoMeshToOpen();
TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh);
TPZCompMesh* CreateCMeshCohe(TPZGeoMesh *gmesh,REAL val);
TPZCompMesh* CreateCMeshToOpen(TPZGeoMesh *gmesh, REAL val);

TPZGeoMesh* CreateGeomeshContDisc();
TPZCompMesh* CreateCmeshContDisc(TPZGeoMesh *gmesh);

void SolveLinearElasticity(TPZCompMesh *cmesh);
void SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an);
void GetSolAtLeft(TPZCompMesh	*cmesh);
void UpdateBcValue(TPZCompMesh *cmesh, REAL val);
  

#endif
