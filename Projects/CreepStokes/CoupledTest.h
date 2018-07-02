/*
 *  CoupledTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__CoupledTest__
#define __PZ__CoupledTest__

#include <cmath>
#include <set>

#include <iostream>
#include <fstream>
#include <string>
#include "pzgmesh.h"
#include "pzstack.h"
#include "TPZVTKGeoMesh.h"
#include "pzanalysis.h"
#include "pzbndcond.h"
#include "TPZCouplingDSMaterial.h"
#include "TPZStokesMaterial.h"
#include "TPZDarcyPMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
#include "pzmat1dlin.h"
#include "pzmat2dlin.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZSkylineNSymStructMatrix.h"
#include "pzstepsolver.h"
#include "TPZGeoLinear.h"
#include "tpzgeoelrefpattern.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSSpStructMatrix.h"
#include "pzanalysis.h"


using namespace std;
using namespace pzshape;

class CoupledTest{
private:
    
    int fdim; //Dimensão do problema
    int fmatID; //Materia do elemento volumétrico
    
    int fmatIdS; //Material Stokes
    int fmatIdD; //Material Darcy
    
    //Materiais das condições de contorno (Stokes)
    int fmatSBCbott;
    int fmatSBCtop;
    int fmatSBCleft;
    int fmatSBCright;
    
    //Materiais das condições de contorno (Darcy)
    int fmatDBCbott;
    int fmatDBCtop;
    int fmatDBCleft;
    int fmatDBCright;
    
    //Material do elemento de interface
    int fmatSInterface;
    int fmatDInterface;
    int fmatInterfaceDS;
    
    //Materiais das condições de contorno (elementos de interface - Stokes)
    int fmatIntSBCbott;
    int fmatIntSBCtop;
    int fmatIntSBCleft;
    int fmatIntSBCright;
    
    //Materiais das condições de contorno (elementos de interface - Darcy)
    int fmatIntDBCbott;
    int fmatIntDBCtop;
    int fmatIntDBCleft;
    int fmatIntDBCright;
    
    //Materia de um ponto
    int fmatPoint;
    int fmatPoint2;
    
    //Condições de contorno do problema
    int fdirichlet;
    int fneumann;
    int fpenetration;
    int fpointtype;
    int fdirichletvar;
    
    
    int fquadmat1; //Parte inferior do quadrado
    int fquadmat2; //Parte superior do quadrado
    int fquadmat3; //Material de interface
    
    STATE fviscosity;
    STATE fpermeability;
    STATE ftheta;
    
    
public:
    
    bool fisH1;
    
    
    CoupledTest();
    
    ~CoupledTest();
    
    void Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE permeability, STATE theta);
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy);
    
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    
    /* Malhas computacionais */
    
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
    
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE permeability, STATE theta);
    
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    // static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    
    static void AddInterfaceCoupllingDS(TPZGeoMesh *gmesh,TPZCompMesh *cmesh, int matInterfaceDS ,int matleft, int matright);
    
};


#endif
