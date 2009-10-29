/*
 *  GenerateSquare.h
 *  SubstructEigen
 *
 *  Created by Philippe Devloo on 29/11/08.
 *  Copyright 2008 UNICAMP. All rights reserved.
 *
 */
#include "tconfig.h"
#include "pzgmesh.h"


TPZGeoMesh *Square(TPZGeoMesh &model, TConfig &cf);

TPZGeoMesh *SquareSingular(TPZGeoMesh &model, TConfig &cf);

TPZGeoMesh *SquareQuad(TPZGeoMesh &model, TConfig &cf);

TPZGeoMesh * DefaultMesh(double R);
TPZGeoMesh * CurvedTriangleMesh(TPZGeoMesh &model, TConfig &cf);
TPZGeoMesh * QuadrilateralMesh(TPZGeoMesh &model, TConfig &cf);

TPZGeoMesh * BlendTriangleMesh(double R);
