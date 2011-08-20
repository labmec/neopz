/**
 * @file
 * @brief Contains the functions to create different computational elements (one- two- three-dimensional).
 */
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

/** @brief Creates computational point element */
TPZCompEl *CreatePointEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational linear element */
TPZCompEl *CreateLinearEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational quadrilateral element */
TPZCompEl *CreateQuadEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational triangular element */
TPZCompEl *CreateTriangleEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational cube element */
TPZCompEl *CreateCubeEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational prismal element */
TPZCompEl *CreatePrismEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational pyramidal element */
TPZCompEl *CreatePyramEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);
/** @brief Creates computational tetrahedral element */
TPZCompEl *CreateTetraEl(TPZGeoEl *gel,TPZCompMesh &mesh,int &index);

#endif