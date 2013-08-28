/*
 *  Tools.h
 *  PZ
 *
 *  Created by Denise de Siqueira on 9/5/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

/*
 * @file
 * @author Denise de Siqueira
 * @since 6/9/11.
 */
#ifndef TOOLSHH
#define TOOLSHH

#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <iostream>
#include <cstdlib>
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "tpzcompmeshreferred.h"
#include "pzanalysis.h"
#include "pzcompel.h"
#include "pzcmesh.h"
#include "pzstepsolver.h"
#include "pzmatrix.h"
#include "pzfstrmatrix.h"
#include "pzlog.h"
#include "pzelctemp.h"
#include "pzfstrmatrix.h"
#include "pzgengrid.h"
#include "pzbndcond.h"
#include "pzmaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "TPZRefPattern.h"


//Criacao de Malhas Geometricas

// Malha Triangular

TPZGeoMesh * MalhaGeoTetraedro( const int h, bool hrefine);
//Malha Quadrilateral formada por um elemento e refinamento nao uniforme
TPZGeoMesh * MalhaGeo( const int h, bool hrefine);
//Malha Quadrilateral formada por 4 elemento e refinamento uniforme
TPZGeoMesh * MalhaGeo2(const int h);
//Malha computacional convencional
TPZCompMeshReferred *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder);
//Malha computacional com p adaptatividade
TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder,bool prefine);
// funcao de forcing funciton
void Forcing1(const TPZVec<REAL> &pt, TPZVec<REAL> &disp);

// definir a funcao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux );

//metodo para graficar a sol no mathematica
void SolGraf(TPZCompMesh *malha, std::ofstream &GraficoSol);

//aplicacao de condicoes de contorno
void CC1(const TPZVec<REAL> &pt, TPZVec<REAL> &f);
void CC2(const TPZVec<REAL> &pt, TPZVec<REAL> &f);
void CC3(const TPZVec<REAL> &pt, TPZVec<REAL> &f);
void CC4(const TPZVec<REAL> &pt, TPZVec<REAL> &f);

//Impressao da malha geomentrica
void PrintGMeshVTK(TPZGeoMesh * gmesh, std::ofstream &file);

//resolve o sistema usando LU
void SolveLU ( TPZAnalysis &an );

void ChangeP(TPZCompMesh * cmesh, TPZCompEl * cel, int newP);

#endif