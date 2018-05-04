/*
 *  StokesTest.cpp
 *  PZ
 *
 *  Created by Pablo Carvalho on 28/07/2017.
 *  Copyright 2017 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef __PZ__HStokesTest__
#define __PZ__HStokesTest__

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
#include "TPZStokesMaterial.h"

#include <pzgeoel.h>
#include "pzgeoelbc.h"
#include "pzfmatrix.h"
#include "pzbstrmatrix.h"
#include <TPZGeoElement.h>
#include "TPZVTKGeoMesh.h"
#include "pzbuildmultiphysicsmesh.h"
#include "TPZInterfaceEl.h"
#include "TPZMultiphysicsInterfaceEl.h"
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

class HStokesTest{
private:
    
    int fdim = 2; //Dimensão do problema
    int fmatID = 1; //Materia do elemento volumétrico
    
    int fmatIDFlux = 2; // material for flux only elements
    //Materiais das condições de contorno
    int fmatBCbott = -1;
    int fmatBCtop = -2;
    int fmatBCleft = -3;
    int fmatBCright = -4;
    
    //Material do elemento de interface
    int fmatInterface = 4;
    
    int ftangentVelocity = 10;
    
    //Materiais das condições de contorno (elementos de interface)
    int fmatIntBCbott = -11;
    int fmatIntBCtop = -12;
    int fmatIntBCleft = -13;
    int fmatIntBCright = -14;
    
    //Material de um ponto
    int fmatPoint = -5;
    
    //Condições de contorno do problema
    int fdirichlet = 0;
    int fneumann = 1;
    int fpenetration = 2;
    int fpointtype = 5;
    int fdirichletvar = 4;
    
    
    int fquadmat1 = 1; //Parte inferior do quadrado
    int fquadmat2 = 2; //Parte superior do quadrado
    int fquadmat3 = 3; //Material de interface
    
    STATE fviscosity = 1.;
    STATE fpermeability = 1.;
    STATE ftheta = -1.;
    
public:

    HStokesTest();
    
    ~HStokesTest();
    
    void Run(int Space, int pOrder, int nx, int ny, double hx, double hy, STATE visco, STATE theta, STATE Sigma);
    
    /*  Malhas geometricas */
    TPZGeoMesh *CreateGMesh(int nx, int ny, double hx, double hy);
    
    //   TPZGeoMesh *GMeshDeformed(int dim, bool ftriang, int ndiv);
    
    
    /* Malhas computacionais */
    
    TPZCompEl *CreateInterfaceEl(TPZGeoEl *gel,TPZCompMesh &mesh,int64_t &index);
    
    TPZCompMesh *CMesh_v(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_p(TPZGeoMesh *gmesh, int Space, int pOrder);
    TPZCompMesh *CMesh_m(TPZGeoMesh *gmesh, int Space, int pOrder, STATE visco, STATE theta, STATE Sigma);
    
    
    //solucao exata
    static void Sol_exact(const TPZVec<REAL> &x, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol);
    
    //lado direito da equacao
    static void F_source(const TPZVec<REAL> &x, TPZVec<STATE> &f, TPZFMatrix<STATE>& gradu);
    
    // static void AddMultiphysicsInterfaces(TPZCompMesh &cmesh, int matfrom, int mattarget);
    void AddMultiphysicsInterfaces(TPZCompMesh &cmesh);
    
    
    
    
};


#endif 
