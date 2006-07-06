//$Id: main.cpp,v 1.1 2006-07-06 15:48:25 tiago Exp $

/**
 * Transient transport equation. Validation test.
 * Jun 06, 2006
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
#include "pztransientanalysis.h"
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
#include "TPZTimer.h"
using namespace std;

int main(){

  TPZCompElDisc::SetOrthogonalFunction(pzshape::TPZShapeDisc::Legendre);

  int p, h;
  cout << "Please, enter h and p parameters: \n";
  cin >> h >> p;
  cout << "\np= " << p << " - h= " << h << "\n";
  TPZCompEl::gOrder = p;
  
  TPZCompMesh * cmesh = CreateMesh(h);    
  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZTransientAnalysis an(cmesh);
  an.TimeStep() = 300.*0.05;
  an.SetNewtonConvergence(1, 1e-14);
  an.SetConvergence(10, 1e-10);
  TPZBandStructMatrix full(cmesh); //  TPZFStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect( ELU );
  an.SetSolver(step);

  cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
   
  char filename[20];
  sprintf(filename,"erro.dat");
  char filedx[20];
  sprintf(filedx,"sol.dx");
  
  an.Run();
   
/**** Aqui faz DX ****/
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);
  an.PostProcess(2);

  an.SetExact(ExactSolution); 
  TPZVec<REAL> pos;
  ofstream out(filename);
  an.PostProcess(pos,out);
  out << "\nNumero de equacoes: " << an.Solution().Rows() << endl;  

  delete cmesh;
  delete gmesh;

  return 0;
}

