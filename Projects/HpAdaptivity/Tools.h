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
#include <pz_config.h>
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
#include "TPZMaterial.h"
#include "tpzquadrilateral.h"
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "pzlog.h"
#include <cmath>
#include "TPZRefPattern.h"


//Criacao de Malhas Geometricas

// Malha Triangular

TPZGeoMesh * MalhaGeoT( const int h, bool hrefine);
//Malha Quadrilateral formada por um elemento e refinamento nao uniforme
TPZGeoMesh * MalhaGeo( const int h, bool hrefine);
//Malha Quadrilateral formada por 4 elemento e refinamento uniforme
TPZGeoMesh * MalhaGeo2(const int h);
TPZGeoMesh *GMesh(bool ftriang, REAL Lx, REAL Ly);
//Malha computacional convencional
TPZCompMeshReferred *CreateCompMesh2d(TPZGeoMesh &gmesh,int porder);
//Malha computacional com p adaptatividade
TPZCompMesh *CompMeshPAdap(TPZGeoMesh &gmesh,int porder,bool prefine);
void SetDifferentOrderP(TPZCompMesh *comp,int porder);
void SetDifferentOrderPMesh4Elem(TPZCompMesh *comp,int porder);
void RefineGeoElements(int dim,TPZGeoMesh *gmesh,TPZVec<REAL> &point,REAL r,REAL &distance,bool &isdefined);
void RefiningNearCircunference(int dim,TPZGeoMesh *gmesh,int nref,int ntyperefs);
void GetPointsOnCircunference(int npoints,TPZVec<REAL> &center,REAL radius,TPZVec<TPZManVector<REAL> > &Points);
void UniformRefine(TPZGeoMesh* gmesh, int nDiv);
void NoUniformRefine(TPZGeoMesh* gmesh, int nDiv);

void RegularizeMesh(TPZGeoMesh *gmesh);
// funcao de forcing funciton
void Forcing1(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingTang(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);

// definir a funcao exata
void SolExata(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux );
void SolExata2(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux );
void SolExata3(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux );
void SolArcTan(const TPZVec<REAL> &pt, TPZVec<STATE> &p, TPZFMatrix<STATE> &flux);


//metodo para graficar a sol no mathematica
void SolGraf(TPZCompMesh *malha, std::ofstream &GraficoSol);
//cria malhas do tipo Grid
TPZGeoMesh * GeoMeshGrid( int h);

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

void ErrorHDiv(TPZCompMesh *hdivmesh, std::ostream &out);
void ErrorL2(TPZCompMesh *l2mesh, std::ostream &out);

//testes do juan
void DirichletEsquerda(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAcima(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void NeumannAbaixo(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void ForcingMista(const TPZVec<REAL> &pt, TPZVec<STATE> &disp);
void SolExataMista(const TPZVec<REAL> &pt, TPZVec<STATE> &solp, TPZFMatrix<STATE> &flux);
TPZCompMesh *CompMeshPAdapJuan(TPZGeoMesh &gmesh,int porder,bool prefine);

#endif
