//$Id: main.cpp,v 1.7 2006-07-06 15:50:34 tiago Exp $

/**
 * Percolation of water from the fracture into the porous media.
 * A property like temperature is convected in the water flow.
 * January 10, 2006
 */

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

#include <time.h>
#include <stdio.h>
#include "meshes.h"
#include "pzcoupledtransportdarcy.h"
#include "pztransientanalysis.h"
#include "TPZTimer.h"

using namespace pzgeom;
using namespace pzshape;
using namespace pzrefine;
using namespace std;
 
/** In this main procedure the coupling of equations is performed by the material TPZCoupledTransportDarcy
 * which behaves like Darcy material in a first moment and then behaves like a transport material using the previous
 * result of Darcy simulation as initial data.
 */
int mainTPZCoupledTransportDarcy(){   

  int h, p;
  cout << "h p\n";
  cin >> h;
  cin >> p;

  TPZTimer WholeTime;
  WholeTime.start();

  TPZMaterial::gBigNumber= 1.e12;
  cout << "h = " << h << " - p = " << p << endl;
  TPZCompMesh * cmesh = /*CheckBetaNonConstant*/CreateSimpleMeshWithExactSolution2(h,p);
  std::cout << "Numero de elementos = " << cmesh->ElementVec().NElements() << std::endl;
  std::cout << "Numero de equacoes  = " << cmesh->NEquations() << std::endl;
  TPZGeoMesh *gmesh = cmesh->Reference();
 
  TPZAnalysis an(cmesh);

//#define ITER_SOLVER
#ifdef ITER_SOLVER
  cout << "ITER_SOLVER" << endl;  
  /*TPZFStructMatrix*/ /*TPZSpStructMatrix*/ TPZBandStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);  
  REAL tol = 1.e-11;
  /*TPZStepSolver*/ TPZCopySolve precond( full.Create() );
  REAL precondtol = 1.E-20;
//  precond.SetSSOR( 1, 1., precondtol, 0);     
//  precond.SetJacobi(1, precondtol, 0);
  step.ShareMatrix( precond );  
  step.SetGMRES( 20000, 30, precond, 1.e-16, 0 ); 
  cout << "SSOR/JACOBI PRECOND" << endl;
  an.SetSolver(step);
#endif  

#define DIRECT
#ifdef DIRECT  
  /*TPZParFrontStructMatrix*/TPZFrontStructMatrix <TPZFrontNonSym>/*TPZFStructMatrix*//*TPZBandStructMatrix*/ full(cmesh);
  full.Create();
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
   step.SetDirect(ELU);
  an.SetSolver(step);
#endif

  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;

  TPZCoupledTransportDarcy::SetCurrentMaterial(0);
  std::cout << "\nCalling an.Run() for FirstEq\n";
  an.Assemble();
  for(int ii = 0; ii < 0*3; ii++)
  {
    TPZMatrix * fmat = step.Matrix();
    ofstream matrizfile("matriz.txt");
    fmat->Print("", matrizfile,EMathematicaInput);
    matrizfile.flush();
    int n = fmat->Rows();
    TPZFMatrix cp(n,n);
    for(int i = 0; i < n; i++) for(int j = 0; j < n; j++) cp(i,j) = fmat->Get(i,j);
    cout << " CN[" << ii << "] = " << cp.ConditionNumber(ii, 20000, 1e-10) << endl;
  }  
  an.Solve();
  an.SetExact(ExactSol_p); //  an.SetExact(SolExata);
  TPZVec<REAL> pos;
  an.PostProcess(pos,std::cout);
  std::cout << "Problem solved\n";
  
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "1stEq_";
  filedx << "Solution.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(4);  
  }      

  TPZCoupledTransportDarcy::SetCurrentMaterial(1);
  std::cout << "\nCalling an.Run() for SecondEq\n";
  an.Solution().Zero();
  an.Assemble();
  an.Solution().Zero();
  an.Solve();
  an.SetExact(ExactSol_u2);
  an.PostProcess(pos,std::cout);  
  std::cout << "Problem solved\n";  
   
  {
  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(0);
  scalnames[0] = "Solution";
//  vecnames[0] = "Derivative";
  std::stringstream filedx;
  filedx << "2ndEq_";
  filedx << "Solution.dx";
  an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
  an.PostProcess(4);  
  }    
    
  WholeTime.stop();
  cout.flush();
  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;
  cout << WholeTime << endl;
  cout << "Banda = " << cmesh->BandWidth() << endl;
  
  delete cmesh;
  delete gmesh;
  return 0; 
}

/** In this main procedure the coupling of equations is performed by a new computational element type: TPZReferredCompEl
 * which can handle solutions between different computational meshes.
 */
#include "meshesReferredCompEl.h"
#include "pzdxmesh.h"
#include "TPZInterfaceEl.h"
int mainTesteConvectivoPuro(){   

  int h, p;
  cout << "h p\n";
//   cin >> h;
//   cin >> p;
  h = 0; p = 2;

  TPZInterfaceElement::SetCalcStiffContDisc();

  TPZTimer WholeTime;
  WholeTime.start();

  TPZMaterial::gBigNumber= 1.e12;
  cout << "h = " << h << " - p = " << p << endl;
  
  TPZCompMesh * firstmesh = /*TesteConvectivoPuro(h, 0, p);*/ TesteConvectivoPuro2(h, p);
  TPZGeoMesh * gmesh = firstmesh->Reference();
  firstmesh->LoadReferences();
  
  if (1){
    ofstream mesh1file("cmesh.txt");
    if (firstmesh) firstmesh->Print(mesh1file);
    ofstream gmeshfile("gmesh.txt");
    if (gmesh) gmesh->Print(gmeshfile);
  }
  
  TPZTransientAnalysis an(firstmesh);
  an.SetConvergence( 50, 1e-10 );
  an.SetNewtonConvergence( 3, 1e-10 );
  an.TimeStep() = 0.1;
  an.SetInitialSolutionAsZero();

  /*TPZParFrontStructMatrix<TPZFrontNonSym>*/TPZFStructMatrix full(firstmesh);
  full.Create();
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);

  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;

  std::cout << "\nCalling an.Run() for FirstEq\n";
  an.Run();
  std::cout << "Problem solved\n";
  
    TPZVec<char *> scalnames(1);
    TPZVec<char *> vecnames(0);
    scalnames[0] = "Solution";
    //vecnames[0] = "Derivative";
    std::stringstream filedx;
    filedx << "convecpuro.dx";
    an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
    an.PostProcess(2);  
    an.CloseGraphMesh();

  
  delete firstmesh;
  delete gmesh;
  return 0; 
}

int main(){   

  int h, p;
  cout << "h p\n";
//   cin >> h;
//   cin >> p;
  h = 4; p = 2;

  TPZInterfaceElement::SetCalcStiffContDisc();

  TPZTimer WholeTime;
  WholeTime.start();

  TPZMaterial::gBigNumber= 1.e12;
  cout << "h = " << h << " - p = " << p << endl;
  
  TPZVec< TPZCompMesh * > CompMeshes(2, NULL);
  
//   CreateMesh_ComoPhilippeQuer_Adimensional_Sem_Simetria(CompMeshes, h, 10, p);
  CreateSimpleMeshesWithExactSolutionToReferredCompEl(CompMeshes, h, p);
  TPZCompMesh * firstmesh  = CompMeshes[0];
  TPZCompMesh * secondmesh = CompMeshes[1];
 
  std::cout << "Numero de elementos = " << CompMeshes[0]->ElementVec().NElements() << std::endl;
  std::cout << "Numero de equacoes  = " << CompMeshes[0]->NEquations() << std::endl;
  TPZGeoMesh *gmesh = CompMeshes[0]->Reference();
   
  firstmesh->LoadReferences();
  
  if (0){
    ofstream mesh1file("cmesh1.txt");
    if (firstmesh) firstmesh->Print(mesh1file);
    ofstream mesh2file("cmesh2.txt");
    if (secondmesh) secondmesh->Print(mesh2file);
    ofstream gmeshfile("gmesh.txt");
    if (gmesh) gmesh->Print(gmeshfile);
  }
  
  TPZAnalysis an(firstmesh);

  TPZParFrontStructMatrix /*TPZFrontStructMatrix*/ <TPZFrontNonSym>/*TPZFStructMatrix*//*TPZBandStructMatrix*/ full(firstmesh);
  full.Create();
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect(ELU);
  an.SetSolver(step);

  std::cout << "Numero de equacoes = " << an.Mesh()->NEquations() << std::endl;

  std::cout << "\nCalling an.Run() for FirstEq\n";
  TPZTimer firstRun;
  firstRun.start();
  an.Run();
  firstRun.stop();
  cout << "FirstRun = " << firstRun << endl;
  an.SetExact(ExactSol_p);
  TPZVec<REAL> pos;
  an.PostProcess(pos,std::cout);
  std::cout << "Problem solved\n";
  
  {
    TPZVec<char *> scalnames(1);
    TPZVec<char *> vecnames(1);
    scalnames[0] = "Solution";
    vecnames[0] = "Derivate";
    std::stringstream filedx;
    filedx << "1stEq_";
    filedx << "Solution.dx";
    an.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
    an.PostProcess(0);  
    an.CloseGraphMesh();
  }
    
  firstmesh->LoadReferences();
  int truefalse = false;
//   cin >> truefalse;
  TPZTransientAnalysis an2(secondmesh, truefalse);
  an2.SetConvergence( 6/*0*/, 1e-10 );
  an2.SetNewtonConvergence( 3, 1e-10 );
  an2.TimeStep() = 50.0;
  an2.SetInitialSolutionAsZero();

  std::cout << "\nCalling an.Run() for SecondEq\n";
  std::cout.flush();

  /*TPZParFrontStructMatrix*/ TPZFrontStructMatrix <TPZFrontNonSym>/*TPZFStructMatrix*//*TPZBandStructMatrix*/ full2(secondmesh);
  full2.Create();
  an2.SetStructuralMatrix(full2);
  TPZStepSolver step2;
  step2.SetDirect(ELU);
  an2.SetSolver(step2);

  an2.Solution().Zero();
  TPZTimer SecondRun;
  SecondRun.start();
  an2.Run();
  SecondRun.stop();
  an2.SetExact(ExactSol_u);
  an2.PostProcess(pos,std::cout);  
  std::cout << "Problem solved\n";  
   
  {
    TPZVec<char *> scalnames(1);
    TPZVec<char *> vecnames(0);
    scalnames[0] = "Solution";
  //  vecnames[0] = "Derivative";
    std::stringstream filedx;
    filedx << "2ndEq_";
    filedx << "Solution.dx";
    an2.DefineGraphMesh(2,scalnames,vecnames,&(filedx.str()[0]));
    an2.PostProcess(2);  
    an2.CloseGraphMesh();
  }    
    
  WholeTime.stop();
  cout.flush();
  std::cout << "Numero de equacoes = " << an2.Mesh()->NEquations() << std::endl;
  cout << WholeTime << endl;
  cout << "Banda = " << firstmesh->BandWidth() << endl;
  cout << "SecondRun time = " << SecondRun << endl;
  
  delete firstmesh;
  delete secondmesh;
  delete gmesh;
  return 0; 
}
