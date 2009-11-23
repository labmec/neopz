#ifndef CREATECONTINUOUSHPP
#define CREATECONTINUOUSHPP
/*
 *  pzcreatecontinuous.h
 *  NeoPZ
 *
 *  Created by Philippe Devloo on 20/11/09.
 *  Copyright 2009 UNICAMP. All rights reserved.
 *
 */
class TPZCompEl;
class TPZCompMesh;
class TPZGeoEl;

TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif