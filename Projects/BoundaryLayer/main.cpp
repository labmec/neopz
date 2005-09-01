//$Id: main.cpp,v 1.4 2005-09-01 19:06:48 tiago Exp $

/**
 * Galerkin descontinuo: problema de camada limite
 * 03/02/2005
 */


#include "meshes.h"
#include "pzvec.h"

#include "pzcmesh.h"

#include "pzdebug.h"
#include "pzcheckgeom.h"
//#include "pzerror.h"

#include "pzgeoel.h"
#include "pzgnode.h"
#include "pzgeoelside.h"

//#include "pzintel.h"
#include "pzcompel.h"
#include "TPZCompElDisc.h"

#include "pzmatrix.h"

#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzbdstrmatrix.h"
#include "pzblockdiag.h"
#include "TPZSpStructMatrix.h"
#include "TPZCopySolve.h"
#include "TPZStackEqnStorage.h"


#include "pzbstrmatrix.h"
#include "pzstepsolver.h"    		
#include "pzonedref.h"

#include "pzadmchunk.h"


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include <time.h>
#include <stdio.h>

#include "gmres.h"

using namespace std;

int main(){


int nmaxeq = 5000;//numero maximo de equacoes do problema
for(int pp = 1; pp < 9; pp++){
 for(int hh = 0; hh < 8; hh++){

//  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int p, h;
  p = pp;
  h = hh;
  TPZCompEl::gOrder = p;
  TPZCompElDisc::gDegree = p;

  TPZCompMesh *cmesh;
//  cmesh = DiscontinuousOnBoundaryLayer(h); 
  cmesh = CreateMeshContDisc(h); 
//  cmesh = CreateMesh(h);
    
  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZAnalysis an(cmesh);

#define direct
#ifdef direct

// Melhor caso para GEM e elementos finitos. Melhor metodo direto
    TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
    //TPZFStructMatrix full(cmesh);
    an.SetStructuralMatrix(full);
    TPZStepSolver step;
    step.SetDirect( ELU /*ECholesky*/ /*ELDLt*/ );
    an.SetSolver(step);

#endif

//#define iter // ACHO QUE A MATRIZ ESPARSA TEM DEFEITO
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix */TPZSpStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);
  
  REAL tol = 1.e-14;

// Melhor caso para Baumann
//      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , true);
//      step.SetGMRES( 200000, 20, *precond, tol, 0 ); 
//      delete precond;

// Sem pre-condicionador 
     TPZCopySolve precond( full.Create() );
     step.ShareMatrix( precond );  
     step.SetGMRES( 2000, 20, precond, tol, 0 ); 

     an.SetSolver(step);
#endif

  cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
   if (an.Mesh()->NEquations() > nmaxeq) { 
     cout << "skipping simulation...\n" << endl;
     delete cmesh;
     delete gmesh;
     break;
   }
  char filename[20];
//  sprintf(filename,"baumann_p%d_h%d.dat",p,h);
//  sprintf(filename,"ef_p%d_h%d.dat",p,h);
  sprintf(filename,"1cont_p%d_h%d.dat",p,h);
  char filedx[20];
//  sprintf(filedx,"baumann_p%d_h%d.dx",p,h);
//  sprintf(filedx,"ef_p%d_h%d.dx",p,h);
  sprintf(filedx,"1cont_p%d_h%d.dx",p,h);

  
  
/*{  
  TPZFMatrix fillin;
  cmesh->ComputeFillIn(50,fillin);
  //fillin.Print("Fillin of the computable mesh");
  VisualMatrix(fillin , filedx);
}*/

  ofstream cmeshout("cmesh.txt");
  cmesh->Print(cmeshout);
  ofstream gmeshout("gmesh.txt");
  gmesh->Print(gmeshout);

  an.Run();

/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);
  an.PostProcess(2);

//   if (gMeshType == 2)
//     an.PostProcess(0);


  an.SetExact(ExactSolution);
  TPZVec<REAL> pos;
  ofstream out(filename);
  an.PostProcess(pos,out);
  out << "\nNumero de equacoes: " << an.Solution().Rows() << endl;  
  
  delete cmesh;
  delete gmesh;
}
}

  return 0;
}

