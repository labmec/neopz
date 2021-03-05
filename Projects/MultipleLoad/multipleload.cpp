
//criando objeto para a classe material
//REAL E = 200000.0, G = 173913.0, v = 0.15;
//REAL E = 10000.0, G = 4347.83, v = 0.15;//teste
//REAL E = 15.133e6, G = 2.59398e7, v = 0.471;//Pinus Pinaster ait.
//REAL E = 69.0e6, G = 5.14378e6, v = 0.330;//Aluminio
//REAL eppx=E, eppy=E, eppz=E, vxy=v, vyz=v,vzx=v, gxy=G, gyz=G, gzx=G;



#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <string.h>
#include <math.h>
#include <stdio.h>

#include "TPZGenGrid2D.h"
#include "pzelasmat.h"
#include "pzbndcond.h"
#include "TPZVTKGeoMesh.h"

#include "pzcmesh.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"
#include "pzstepsolver.h"

#include "pzlog.h"

int gDebug;

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multipleload.main"));
#endif

/**
 * -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <- -> PROGRAMA PRINCIPAL <-
 */

TPZAutoPointer<TPZGeoMesh> GenerateGMesh();
TPZAutoPointer<TPZCompMesh> GenerateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh);

int main()
{
    InitializePZLOG();

    TPZAutoPointer<TPZGeoMesh> gmesh = GenerateGMesh();
    TPZAutoPointer<TPZCompMesh> cmesh = GenerateCMesh(gmesh);
    
    TPZAnalysis an(cmesh);
    TPZSkylineStructMatrix skylstr(cmesh);
    an.SetStructuralMatrix(skylstr);
    TPZStepSolver<STATE> step;
    step.SetDirect(ECholesky);
    an.SetSolver(step);
    an.Run();
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Solution().Print("Solution",sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    std::string post1("load.0.vtk");
    TPZStack<std::string> scalnames,vecnames;
    vecnames.Push("displacement");
    scalnames.Push("SigmaY");
    an.DefineGraphMesh(2, scalnames, vecnames, post1);
    an.PostProcess(2);
    TPZElasticityMaterial *mat = dynamic_cast<TPZElasticityMaterial *> (cmesh->MaterialVec()[1]);
    mat->SetPostProcessIndex(1);
    an.PostProcess(2);
    
    
    return 0;
}

TPZAutoPointer<TPZCompMesh> GenerateCMesh(TPZAutoPointer<TPZGeoMesh> gmesh)
{
    TPZAutoPointer<TPZCompMesh> cmesh = new TPZCompMesh(gmesh);
    REAL Elas(100.),Poiss(0.2),fx(0.),fy(0.);
    TPZElasticityMaterial *elasmat = new TPZElasticityMaterial(1,Elas,Poiss,fx,fy);
    cmesh->InsertMaterialObject(elasmat);
    TPZFMatrix<STATE> val1(2,2,0.),val2(2,1,0.);
    TPZVec<TPZFMatrix<STATE> > val2vec(2);
    TPZBndCond *bc2 = new TPZBndCond(elasmat,-2,1,val1,val2);
    val2vec[1] = val2;
    val2(1,0) = 1.;
    val2vec[0] = val2;
    bc2->SetLoadCases(val2vec);
    cmesh->InsertMaterialObject(bc2);
    TPZBndCond *bc3 = new TPZBndCond(elasmat,-3,1,val1,val2);
    val2.Zero();
    val2vec[0] = val2;
    val2(1,0) = 1.;
    val2vec[1] = val2;
    bc3->SetLoadCases(val2vec);
    cmesh->InsertMaterialObject(bc3);
    val1(1,1) = 1.e8;
    val2.Zero();
    TPZBndCond *bc4 = new TPZBndCond(elasmat,-4,2,val1,val2);
    cmesh->InsertMaterialObject(bc4);
    val1.Zero();
    val1(0,0) = 1.e8;
    TPZBndCond *bc5 = new TPZBndCond(elasmat, -5, 2, val1, val2);
    cmesh->InsertMaterialObject(bc5);
    
    cmesh->AutoBuild();
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        cmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    
    return cmesh;
}

TPZAutoPointer<TPZGeoMesh> GenerateGMesh()
{
    //criando objetos das classes computacional e geometrica
    REAL dx = 1.0;
    REAL dy = 1.0;
    int nelx = 10;
    int nely = 10;
 
    TPZManVector<int,3> nx(2,1);
    TPZManVector<REAL,3> x0(3,0.),x1(3,0.);
    x1[0] = dx;
    x1[1] = dy;
    nx[0] = nelx;
    nx[1] = nely;
    TPZGenGrid2D gengrid(nx,x0,x1);
    TPZAutoPointer<TPZGeoMesh> gmesh = new TPZGeoMesh;
    gengrid.Read(gmesh);
    x1[1] = 0.;
    x1[0] = 0.2;
    gengrid.SetBC(gmesh, x0, x1, -2);
    x0 = x1;
    x1[1] = 0.;
    x1[0] = 0.5;
    gengrid.SetBC(gmesh, x0, x1, -3);
    x0 = x1;
    x1[0] = 1.;
    gengrid.SetBC(gmesh, x0, x1, -4);
    x0[0] = 0.;
    x0[1] = 1.;
    x1[0] = 0.;
    x1[1] = 0.;
    gengrid.SetBC(gmesh, x0, x1, -5);
    
#ifdef LOG4CXX
    {
        std::stringstream sout;
        gmesh->Print(sout);
        LOGPZ_DEBUG(logger, sout.str())
    }
#endif
    
    std::ofstream out("Geometry.vtk");
    TPZVTKGeoMesh::PrintGMeshVTK(gmesh.operator->(), out, true);
    return gmesh;
}
