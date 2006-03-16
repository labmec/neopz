//$Id: main.cpp,v 1.2 2006-03-16 13:26:01 tiago Exp $

/**
 * Validation test of TPZElasticity3D material
 * 24/11/2003
 */

#include "performance.h"
#include "StrMatrix.h"

#include <sstream>
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


#include "pzbndcond.h"
#include "pzpoisson3d.h"

#include "pzvisualmatrix.h"

#include "TPZGeoElement.h"
#include "pzgeoel.h"

#include "TPZGeoCube.h"
#include "pzshapecube.h"
#include "TPZRefCube.h"
#include "pzshapelinear.h"
#include "TPZGeoLinear.h"
#include "TPZRefLinear.h"
#include "pzrefquad.h"
#include "pzshapequad.h"
#include "pzgeoquad.h"
#include "pzshapetriang.h"
#include "pzreftriangle.h"
#include "pzgeotriangle.h"
#include "pzshapeprism.h"
#include "pzrefprism.h"
#include "pzgeoprism.h"
#include "pzshapetetra.h"
#include "pzreftetrahedra.h"
#include "pzgeotetrahedra.h"
#include "pzshapepiram.h"
#include "pzrefpyram.h"
#include "pzgeopyramid.h"
#include "pzrefpoint.h"
#include "pzgeopoint.h"
#include "pzshapepoint.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"

#include <time.h>
#include <stdio.h>
#include "meshes.h"
#include "tetgenmesh.h"
#include "TPZTimer.h"
#include "pzlog.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;

#include <set>

enum solver {CG = 0, GMRES = 1, CHOLESKY, LU};
enum strmatrix {SPARSE = 0, SKYLINE = 1, SKYLINEPAR, FRONTAL, FULL};
enum precond {NONE = 0, NODECENTERED = 1, ELEMENT = 2, BLOCKJACOBI, SSOR, JACOBI};



int main(){

  InitializePZLOG();

  const int numiter = 2;
  const REAL tol = 1.e-10;
  const int SSOR_numiter = 10;
  const REAL SSOR_tol = 1.e-10;
  
  solver SOLVER;
  strmatrix STRMATRIX;
  precond PRECOND;

  ifstream file("config.txt");
  ofstream out("result.txt");
  //Read problem parameters h and p. h indicates number of uniform refinements and p is the order of interpolation (constant in the mesh)
  int inteiro;
  file >> inteiro;
  const int h = inteiro;
  file >> inteiro;
  const int p = inteiro;
  
  //Init problem
  TPZCompMesh * cmesh = /*Small(h, p);*/CreateComputeMeshFromTetGen/*BarraTracionadaGirada*//*BarraTracionadaNeumann*/(h, p);
  cout << "h = " << h << ", p = " << p << "\nNEquations = " << cmesh->NEquations() << endl;
  out  << "h = " << h << ", p = " << p << "\nNEquations = " << cmesh->NEquations() << endl;
  TPZGeoMesh * gmesh = cmesh->Reference(); 
  TPZAnalysis an(cmesh);    
  
  //Read solver
  string buf;
  file >> buf;
  if (buf == "CG"){ 
    SOLVER  = CG;
    cout << "CG\n";
    out << "CG\n";
  }
  if (buf == "GMRES"){
    SOLVER = GMRES;
    cout << "GMRES\n";
    out << "GMRES\n";
  }
  if (buf == "CHOLESKY"){
    SOLVER = CHOLESKY;
    cout << "CHOLESKY\n";
    out << "CHOLESKY\n";
  }
  if (buf == "LU"){
    SOLVER = LU;
    cout << "LU\n";
    out << "LU\n";
  }
 
  //Read strmatrix enum strmatrix {SPARSE = 0, SKYLINE = 1, SKYLINEPAR, FRONTAL, FULL}
  TPZStructMatrix * StrMatrix = 0;
  file >> buf;
  if (buf == "SPARSE"){
    STRMATRIX = SPARSE;
    StrMatrix = new TPZSpStructMatrix(cmesh);
    cout << "SPARSE\n";
    out << "SPARSE\n";
  }
  if (buf == "SKYLINE"){
    STRMATRIX = SKYLINE;
    StrMatrix = new TPZSkylineStructMatrix(cmesh);
    cout << "SKYLINE\n";
    out << "SKYLINE\n";
  }
  if (buf == "SKYLINEPAR"){
    STRMATRIX = SKYLINEPAR;
    StrMatrix = new TPZParSkylineStructMatrix(cmesh);
    cout << "SKYLINEPAR\n";
    out << "SKYLINEPAR\n";
  }    
  if (buf == "FRONTAL"){
    STRMATRIX = FRONTAL;
    StrMatrix = new TPZParFrontStructMatrix<TPZFrontSym>(cmesh);
    cout << "FRONTAL\n";
    out << "FRONTAL\n";
  }    
  if (buf == "FULL"){
    STRMATRIX = FULL;
    StrMatrix = new TPZFStructMatrix(cmesh);
    cout << "FULL\n";
    out << "FULL\n";
  }    
  cout.flush();
  
  an.SetStructuralMatrix(*StrMatrix);
  TPZStepSolver step(StrMatrix->Create());
  an.SetSolver(step);


  //ITERATIVE SOLVERS
  if (SOLVER == CG || SOLVER == GMRES){
//    an.Solver().Matrix()->Print("", cout, EMathematicaInput);
    while (file){
      //PRECOND enum precond {NONE = 0, NODECENTERED = 1, ELEMENT = 2, BLOCKJACOBI, SSOR, JACOBI};
      file >> buf;
      if (!file) break;
      TPZMatrixSolver * preconditioner = 0;
      if (buf == "NONE"){
        PRECOND = NONE;
        cout << "Precond = NONE\n";
        out << "Precond = NONE\n";
        preconditioner = new TPZCopySolve( NULL );
      }
      if (buf == "NODECENTERED"){
        PRECOND = NODECENTERED;
        cout << "Precond = NODECENTERED\t";
        out << "Precond = NODECENTERED\t";
        file >> buf;
        bool overlap = false;
        if (buf == "FALSE"){
          out << "FALSE\n";
          cout << "FALSE\n";
          overlap = false;
        }
        if (buf == "TRUE"){ 
          out << "TRUE\n";
          cout << "TRUE\n";
          overlap = true;        
        }
        std::cout << "\nBuilding PreConditioner";
        preconditioner = an.BuildPreconditioner(TPZAnalysis::ENodeCentered , overlap);
        std::cout << "\ndone\n";
      }
      if (buf == "ELEMENT"){
        PRECOND = ELEMENT;
        cout << "Precond = ELEMENT\t";
        out << "Precond = ELEMENT\t";
        file >> buf;
        bool overlap = false;
        if (buf == "FALSE"){
          out << "FALSE\n";
          cout << "FALSE\n";
          overlap = false;
        }
        if (buf == "TRUE"){ 
          out << "TRUE\n";
          cout << "TRUE\n";
          overlap = true;        
        }
        std::cout << "\nBuilding PreConditioner";
        preconditioner = an.BuildPreconditioner(TPZAnalysis::EElement , overlap);
        std::cout << "\ndone\n";
      }      
      if (buf == "BLOCKJACOBI"){
        PRECOND = BLOCKJACOBI;
        cout << "Precond = BLOCKJACOBI\t";
        out << "Precond = BLOCKJACOBI\t";
        file >> buf;
        bool overlap = false;
        if (buf == "FALSE"){
          out << "FALSE\n";
          cout << "FALSE\n";
          overlap = false;
        }
        if (buf == "TRUE"){ 
          out << "TRUE\n";
          cout << "TRUE\n";
          overlap = true;        
        }
        std::cout << "\nBuilding PreConditioner";
        preconditioner = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , overlap);
        std::cout << "\ndone\n";
      }   
      if (buf == "SSOR"){
        PRECOND = SSOR;
        cout << "Precond = SSOR\n";
        out << "Precond = SSOR\n";
        preconditioner = new TPZStepSolver( StrMatrix->Create() );
        dynamic_cast<TPZStepSolver*>(preconditioner)->SetSSOR(SSOR_numiter, 1., SSOR_tol, 0);
      }
      if (buf == "JACOBI"){
        PRECOND = JACOBI;
        cout << "Precond = JACOBI\n";
        out << "Precond = JACOBI\n";
        preconditioner = new TPZStepSolver( StrMatrix->Create() );
        dynamic_cast<TPZStepSolver*>(preconditioner)->SetJacobi(SSOR_numiter, SSOR_tol, 0);
      }
      
      if (SOLVER == CG){
        step.SetCG(numiter, *preconditioner, tol, 0);
      }
      if (SOLVER == GMRES){
        step.SetGMRES( numiter, 20, *preconditioner, tol, 0);
      }
      
      an.SetSolver(step);  
      
      an.Assemble();  
      
      TPZPerformance perf;
      cout << "\nSTART ...\n";
      cout.flush();
              
              an.Solution().Zero();              
              perf.Start();   
              an.Solve();
              perf.Stop();
              
      cout << "\nFINISHED!\n";
        cout.flush();
        out.flush();
      out  << "NITER: " << an.Solver().fNITER << endl;
      cout << "NITER: " << an.Solver().fNITER << endl;
      perf.Print(cout);
      perf.Print(out);
      
      if (preconditioner) delete preconditioner;
      out << "\n\n";
    }//while    
  }//ITERATIVE SOLVERS
  
  if (SOLVER == CHOLESKY || SOLVER == LU){
  
    if (SOLVER == CHOLESKY) step.SetDirect(ECholesky);
    if (SOLVER == LU) step.SetDirect(ELU);
    if (STRMATRIX == FRONTAL){
      cout << "FRONTAL" << endl;
      TPZPerformance perf;
      an.SetSolver(step);  
      perf.Start();
      an.Run();
      perf.Stop();  
      perf.Print(cout);
      perf.Print(out);
    }
    else{
      cout << "DIRECT" << endl;
      an.SetSolver(step);  
      TPZPerformance perf;
      perf.Start();      
      an.Run();
      perf.Stop();  
      perf.Print(cout);
      perf.Print(out);
    }
  }//DIRECT SOLVERS

  
   
  delete cmesh;
  delete gmesh;
  return 0;



}//main




int main22(){
   
//   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Chebyshev;
//  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Jacobi;  
//  TPZShapeLinear::fOrthogonal = TPZShapeLinear::Legendre;  

  cout << "h p" << endl;
  int h, p;
  cin >> h;
  cin >> p;
  cout << "h = " << h << " - p = " << p << endl;  

  TPZCompMesh * cmesh = BarraTracionadaGirada/*BarraTracionadaNeumann*/(h, p);
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

//#define direct  
#ifdef direct
  /*TPZFrontStructMatrix <TPZFrontSym> */TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ECholesky);
  an.SetSolver(step);
#endif
  
#define iter
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix*/ TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-8;

#define Precond
#ifdef Precond 

  #define SSORPRECOND
  #ifdef SSORPRECOND
     TPZStepSolver precond( full.Create() );
     REAL precondtol = 1.E-20;
//     precond.SetSSOR( 1, 1., precondtol, 0);     
     precond.SetJacobi(1, precondtol, 0);
     step.ShareMatrix( precond );  
     step.SetCG( 20000, precond, tol, 0 );
//     step.SetGMRES( 2000, 20, precond, tol, 0);     
     cout << "SSOR/JACOBI PRECOND" << endl;
  #else

      TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EBlockJacobi , false);
      step.SetCG( 2000, *precond, tol, 0 );
//      step.SetGMRES( 2000, 20, *precond, tol, 0);
      delete precond;
  #endif

#else
// Sem pre-condicionador 
     TPZCopySolve precond( full.Create() );
     step.ShareMatrix( precond );  
     step.SetCG( 2000, precond, tol, 0 );
//     step.SetGMRES( 2000, 20, precond, tol, 0);     
//     step.SetSSOR( 18000, 1., tol, 0);
     cout << "SEM PRECOND" << endl;
#endif

     an.SetSolver(step);
#endif  
  
  an.Run();
  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;
  cout << "h = " << h << " - p = " << p << endl;  
  std::cout.flush();
   
  delete cmesh;
  delete gmesh;
  return 0;
}


