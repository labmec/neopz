//
//  TPZCohesiveTests.h
//  PZ
//
//  Created by Agnaldo Farias on 9/5/12.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#ifndef CohesiveTests
#define CohesiveTests


class TPZGeoMesh;
class TPZCompMesh;
class TPZNonLinearAnalysis;

void ElastTest();
void ElastNLTest();
TPZGeoMesh* CreateGeoMesh();
TPZGeoMesh* CreateGeoMeshCohe();
TPZCompMesh* CreateCMesh(TPZGeoMesh *gmesh);
TPZCompMesh* CreateCMeshCohe(TPZGeoMesh *gmesh);
void SolveLinearElasticity(TPZCompMesh *cmesh);
void SolveNLElasticity(TPZCompMesh *cmesh, TPZNonLinearAnalysis &an);
void ElastNLTestWithCohesive();
void GetSolAtLeft(TPZCompMesh	*cmesh);
  

#endif
