//Classes utilitarias
#include "pzvec.h"
#include "pzgmesh.h"
#include "pzgeoel.h"
#include "pzcmesh.h"
#include "pzcompel.h"
#include "pzfmatrix.h"

//Bibliotecas std, math etc
#include <stdlib.h>
#include <iostream>
#include <stdio.h>
#include <time.h>
#include <math.h>

// Classes to computational analysis
#include "TPZCompElDisc.h"
#include "TPZShapeDisc.h"
#include "pzl2projection.h"
#include "tpzautopointer.h"
#include "pzanalysis.h"
#include "TExtFunction.h"
#include "pzlog.h"
#include "pzbstrmatrix.h"

// Enabling namespace
using namespace std;
using namespace pzshape;

/** Acrescenta uma funcao de forma singular */
int main() {
  TPZGeoMesh *gmesh = new TPZGeoMesh;
  REAL co[4][2] = {{-0.5,-0.2},{2.,-0.2},{2.,3.},{-0.5,3.}};
  for(int nod = 0; nod < 4; nod++){
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }
  cout << "\nTriangulo (0) ou Quadrado (1)? ";
  int opcao;
  cin >> opcao;
  int nnodes;
  if (opcao == 0) nnodes = 3;
  if (opcao == 1) nnodes = 4;
  TPZVec<int> nodind(nnodes);
  for(int i = 0; i < nnodes; i++) nodind[i] = i;
  int index;
  TPZGeoEl * geo;
  if (nnodes == 3) geo = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  if (nnodes == 4) geo = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  gmesh->BuildConnectivity();
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  TPZVec<REAL> sol(1,0.);
  TPZAutoPointer<TPZMaterial> material = new TPZL2Projection(1,2,1,sol);
  cmesh->InsertMaterialObject(material);
  cmesh->SetAllCreateFunctionsDiscontinuous();
  TPZCompEl::SetgOrder(1);
  cmesh->SetDefaultOrder(1);
  cmesh->AutoBuild();
  cout << "\ncmesh->NElements() = " << cmesh->NElements() << "\n";
  TPZCompElDisc * cel = dynamic_cast<TPZCompElDisc*>(cmesh->ElementVec()[0]);
  TPZAutoPointer<TPZFunction> extShapes = new TExtFunction();
  cel->SetExternalShapeFunction(extShapes);  
  TPZVec<int> ord(3,10);
	TPZAutoPointer<TPZIntPoints> points = cel->GetIntegrationRule().Clone();
  points->SetOrder(ord);
  const int npoints = cel->GetIntegrationRule().NPoints();
  TPZVec<REAL> qsi(2);
  REAL w;
  TPZFMatrix phi, dphi;
  ofstream saidaX("saidaX.txt");
  ofstream saidaY("saidaY.txt");
  saidaX << "{";
  saidaY << "{";
  for(int i = 0; i < npoints; i++){
    points->Point(i, qsi, w);
    cel->Shape(qsi, phi, dphi);
    TPZVec<REAL> X(3);
    cel->Reference()->X(qsi,X);
    int pos = phi.Rows()-1;
    saidaX << "{" << X[0] << "," << phi(pos,0) << "},";
    saidaY << "{" << X[1] << "," << phi(pos,0) << "},";
  }
  bool dx = true;
  if (dx){
    cmesh->ExpandSolution();
    cmesh->Solution().Zero();
//    int pos = phi.Rows()-1;
    TPZConnect & np = cel->Connect(0);
    int blocknumber = np.SequenceNumber();
    int firsteq = cmesh->Block().Position(blocknumber);
    int ndf = cmesh->Block().Size(blocknumber);
    cmesh->Solution()(firsteq+ndf-1,0) = 1.0;
    
    TPZAnalysis an(cmesh);
    an.Solution() = cmesh->Solution();
    TPZVec<std::string> scalnames(1);
    scalnames[0] = "Solution";
    TPZVec<std::string> vecnames(0);
    char plotfile[] = "singular.dx";
    an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
    an.PostProcess(8);  
  }
  delete cmesh;
  delete gmesh;
  return 0;
}

/** Disco com furo no meio - problema de Poisson */
#include "pzpoisson3d.h"
#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "pzfstrmatrix.h"
#include "pzskylstrmatrix.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZParFrontMatrix.h"
#include "TPZFrontNonSym.h"
#include "pzstepsolver.h"
#include "TDiscoFunction.h"
#include "TPZCopySolve.h"
#include "TPZSpStructMatrix.h"

void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u.Resize(1);
  deriv.Resize(2,1);
  double Xp = x[0];
  double Yp = x[1];
  double r = sqrt( pow(Xp,2) + pow(Yp,2) );
  if(r < 1e-6) r = 1e-6;
  
  /*u[0] = Xp*Xp +Yp*Yp;
  deriv(0,0) = 2.*Xp;
  deriv(1,0) = 2.*Yp;
  return;*/
  
  u[0] = 8.148974294721926 + 2.6704656055626725*log(r);
  deriv(0,0) = (2.6704656055626725*Xp)/(r*r);
  deriv(1,0) = (2.6704656055626725*Yp)/(r*r);
}

void Dirichlet(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
  TPZFNMatrix<2> deriv(2,1);
    TPZManVector<REAL, 3> xcopy(x);
  ExactSolution(xcopy,f,deriv);
}

void LoadFunction(const TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f.Resize(1);
  f[0] = 0.;
}

void UniformRefinement(const int dim, TPZGeoMesh &gmesh, bool AllMatId, const int WhichMatId){
  TPZManVector<TPZGeoEl*> filhos;
  int n = gmesh.NElements();
  for(int i = 0; i < n; i++){
    TPZGeoEl * gel = gmesh.ElementVec()[i];
    if(!gel) continue;
    if(gel->Dimension() != dim) continue;
    if(gel->HasSubElement()) continue;
    if(AllMatId == false){
      if(gel->MaterialId() == WhichMatId){
        gel->Divide(filhos);
      }
    }
    else{
      gel->Divide(filhos);
    }
  }
}

void LocalRefinement(const int dim, const int matid, TPZGeoMesh &gmesh){
  TPZManVector<TPZGeoEl*> filhos;
  int n = gmesh.NElements();
  for(int i = 0; i < n; i++){
    TPZGeoEl * gel = gmesh.ElementVec()[i];
    if(!gel) continue;
    if(gel->Dimension() != dim) continue;
    if(gel->HasSubElement()) continue;
    bool refine = false;
    for(int is = 0; is < gel->NSides(); is++){
      TPZGeoElSide neig = gel->Neighbour(is);
      while(neig.Element() != gel){
        if(neig.Element()->MaterialId() == matid){
          refine = true;
          break;
        }
        neig = neig.Neighbour();
      }
    }
    if (refine) gel->Divide(filhos);
  }
}

void PRefinement(TPZCompMesh &cmesh, const int InitialP){
  int nel = cmesh.NElements();
  for(int iel = 0; iel < nel; iel++){
    TPZCompEl * cel = cmesh.ElementVec()[iel];
    if(!cel) continue;
    TPZInterpolationSpace * sp = dynamic_cast<TPZInterpolationSpace*>(cel);
    if(!sp) continue;
    int neworder = sp->Reference()->Level() + InitialP;
    sp->PRefine( neworder );
  }//for
}//void

int main1(){

  const bool ApenasPolinomial = false;

  InitializePZLOG("log4cxx.cfg");
  TPZShapeDisc::fOrthogonal = TPZShapeDisc::Legendre;

  const int nnodes = 20;
  const double scaleFuro = 1e-1;
  REAL co[nnodes][2] = {{2., 0.}, {1.618033988749895,1.1755705045849463}, {0.6180339887498949, 1.902113032590307}, 
   {-0.6180339887498949, 1.902113032590307}, {-1.618033988749895, 1.1755705045849463}, {-2., 0.},
   {-1.618033988749895, -1.1755705045849463},{-0.6180339887498949, -1.902113032590307}, {0.6180339887498949,-1.902113032590307},
   {1.618033988749895, -1.1755705045849463},
   
    {0.1*scaleFuro, 0.*scaleFuro}, {0.08090169943749476*scaleFuro , 0.058778525229247314*scaleFuro }, {0.030901699437494747*scaleFuro ,   0.09510565162951536*scaleFuro }, 
    {-0.030901699437494747*scaleFuro ,   0.09510565162951536*scaleFuro }, {-0.08090169943749476*scaleFuro ,   0.058778525229247314*scaleFuro }, 
    {-0.1*scaleFuro ,   0.*scaleFuro }, {-0.08090169943749476*scaleFuro , -0.058778525229247314*scaleFuro }, {-0.030901699437494747*scaleFuro , -0.09510565162951536*scaleFuro }, 
    {0.030901699437494747*scaleFuro , -0.09510565162951536*scaleFuro }, {0.08090169943749476*scaleFuro , -0.058778525229247314*scaleFuro }};
   
  const int nelem = 10; 
  int indices[nelem][4] = {{0,1,11,10},{1,2,12,11},{2,3,13,12},{3,4,14,13},{4,5,15,14},{5,6,16,15},{6,7,17,16},{7,8,18,17},{8,9,19,18},{9,0,10,19}};
  const int nelbc = 20;
  int contorno[nelbc][3] = {{0,1,-2},{1,2,-2},{2,3,-2},{3,4,-2},{4,5,-2},{5,6,-2},{6,7,-2},{7,8,-2},{8,9,-2},{9,0,-2},{10,11,-3},{11,12,-3},{12,13,-3},
                            {13,14,-3},{14,15,-3},{15,16,-3},{16,17,-3},{17,18,-3},{18,19,-3},{19,10,-3}};

/** teste - apenas 1 elemento **/
/*  const int nelem = 1; 
  int indices[nelem][4] = {{0,1,11,10}};
  const int nelbc = 2;
  int contorno[nelbc][3] = {{0,1,-2},{10,11,-3}};*/
/** teste - apenas 1 elemento **/


  TPZGeoMesh *gmesh = new TPZGeoMesh();
  for(int nod=0; nod<nnodes; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  const int matid = 1;
  for(int el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(int nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    gmesh->CreateGeoElement(EQuadrilateral,nodind,matid,index);
  }

  for(int el=0; el<nelbc; el++) {
    TPZManVector<int,2> nodind(2);
    for(int nod=0; nod<2; nod++) nodind[nod]=contorno[el][nod];
    int bcmatid = contorno[el][2];
    int index;
    gmesh->CreateGeoElement(EOned,nodind,bcmatid,index);
  }


  gmesh->BuildConnectivity();
  
  int p;
  if(ApenasPolinomial) p = 2;
  else p = 0;
  cout << "npassos = "; int hval; cin >> hval;
  const int h = hval;   cout << "\nh=" << h;
  for(int ir = 0; ir < h; ir++){
//     UniformRefinement(2,*gmesh, true, 8511965); UniformRefinement(1,*gmesh, true, 8511965);
     LocalRefinement(2,-3,*gmesh); UniformRefinement(1,*gmesh, false, -3);
  }
  TPZCompMesh * cmesh = new TPZCompMesh(gmesh);
  
  cmesh->SetDimModel(2);
  
  TPZAutoPointer<TPZMaterial> mat = new TPZMatPoisson3d(matid, 2);
  mat->SetForcingFunction( new TPZDummyFunction(LoadFunction) );
  TPZMatPoisson3d * matcast = dynamic_cast<TPZMatPoisson3d*>(mat.operator->());
  matcast->fPenaltyConstant = 1.;
  matcast->SetSolutionPenalty(); 
//    matcast->SetFluxPenalty();
//   matcast->SetBothPenalty();
//   matcast->SetSymmetric();

  TPZManVector<REAL,2> convdir(2,0.);
  REAL diff = 1.;
  matcast->SetParameters(diff, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,10.);
  TPZAutoPointer<TPZMaterial> bcFora ( mat->CreateBC(mat,-2, 0,val1,val2) );
  val2(0,0) = 2.;
  TPZAutoPointer<TPZMaterial> bcDentro = mat->CreateBC(mat,-3, 0,val1,val2);
  
  bcFora->SetForcingFunction(new TPZDummyFunction(Dirichlet));
  bcDentro->SetForcingFunction(new TPZDummyFunction(Dirichlet));
  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bcFora);
  cmesh->InsertMaterialObject(bcDentro);
  
  cmesh->SetAllCreateFunctionsDiscontinuous();

  TPZCompEl::SetgOrder(p);
  cmesh->SetDefaultOrder(p);
  cmesh->AutoBuild();
  cmesh->ComputeNodElCon();
//   cmesh->AdjustBoundaryElements();

  ofstream malhas("malhas.txt");
  cmesh->Reference()->Print(malhas);
  cmesh->Print(malhas);

  if(ApenasPolinomial){
    PRefinement(*cmesh, 2);
  }
  cmesh->AdjustBoundaryElements();


//   TPZCompElDisc::SetTotalOrderShape(cmesh);
  TPZCompElDisc::SetTensorialShape(cmesh);

  if(!ApenasPolinomial){
    TPZAutoPointer<TPZFunction> ExternalShapes = new TDiscoFunction();
    for(int iel = 0; iel < cmesh->NElements(); iel++){
      TPZCompEl * cel = cmesh->ElementVec()[iel];
      if(!cel) continue;
      TPZCompElDisc * disc = dynamic_cast<TPZCompElDisc*>(cel);
      if(!disc) continue;
      disc->SetExternalShapeFunction(ExternalShapes);
    }
  }//if


  TPZAnalysis an(cmesh);
//   TPZFrontStructMatrix<TPZFrontNonSym> matrix(cmesh);
//   TPZFrontStructMatrix<TPZFrontSym> matrix(cmesh);
//   TPZFStructMatrix matrix(cmesh);
  TPZBandStructMatrix matrix(cmesh);
  an.SetStructuralMatrix(matrix);
  TPZStepSolver step;
  
// #define DIRETO
#ifdef DIRETO  
    step.SetDirect(ELU);//ECholesky); 
#else 
//       TPZCopySolve precond( matrix.Create() );step.ShareMatrix( precond );  
	TPZAutoPointer<TPZGuiInterface> gui;
	TPZFMatrix fakeRhs(cmesh->NEquations(),1);
	TPZBandStructMatrix PrecondMatrix(cmesh); 
	TPZStepSolver precond(PrecondMatrix.CreateAssemble(fakeRhs,gui));
	precond.SetDirect(ELU);
      
      step.SetGMRES( 300000, 160, precond, 1.e-14, 0 );
#endif  
  
  
  an.SetSolver(step);  
  cout << "\nNumero de equacoes = " << cmesh->NEquations() << std::endl;
  cout << "Banda = " << cmesh->BandWidth() << "\n";
  cout.flush();
  an.Assemble();

{ofstream fileKsdvcxcxv("k.nb");
 an.Solver().Matrix()->Print("m=",fileKsdvcxcxv,EMathematicaInput);
 fileKsdvcxcxv << "\nRe[Eigenvalues[m]]\n";
 }
 
  /** checando residuo da sol exata descontinua */
 /* {
            TPZFMatrix mySol = an.Rhs();
            int nshapeperel = TPZShapeDisc::NShapeF(p,2,TPZShapeDisc::EOrdemTotal)+ExternalShapes->NFunctions();
            const int ndiscel = mySol.Rows() / nshapeperel;
            mySol.Zero();
            for(int i = 0; i < ndiscel; i++){
              int pos = (i+1)*nshapeperel - 1; 
              mySol(pos, 0) = 2.6704656055626725;
              pos = (i+1)*nshapeperel - 2;
              mySol(pos, 0) = 8.148974294721926;
            }
            TPZFMatrix myRhs;
            an.Solver().Matrix()->Multiply(mySol, myRhs);
            myRhs -= an.Rhs();
            ofstream residuofile("residuo.nb");
            myRhs.Print("myRhs=",residuofile);
            ofstream rigfile("rigidez.nb");
            an.Solver().Matrix()->Print("rigidez=",rigfile,EMathematicaInput);
            an.Rhs().Print("carga=",rigfile,EMathematicaInput);
            an.Solution() = mySol;
            an.LoadSolution();
  }*/
  /** ate aqui */
  
  
  an.Solve();
  an.LoadSolution();
  
/*  TPZVec<REAL> error(cmesh->NElements(),0.);
  TPZCompElDisc::EvaluateSquareResidual2D(*cmesh, error,false);
  {ofstream resfile("res.txt");
  for(int i = 0; i < error.NElements(); i++) resfile << sqrt(error[i]) << "\n";}
  cmesh->SetElementSolution(0, error);*/

  
      an.SetExact(ExactSolution);
      TPZVec<REAL> pos;
      ofstream errorfile("erro.txt");
      an.PostProcess(pos,errorfile);
      errorfile << "Numero de equacoes = " << cmesh->NEquations() << std::endl;
      errorfile << "Banda = " << cmesh->BandWidth() << "\n";
      errorfile.flush();  
  
//   ofstream solfile("solucao.txt");
//   an.Solution().Print("solucao",solfile);
  
  TPZVec<std::string> scalnames(4);
  scalnames[0] = "Solution";
  scalnames[1] = "p";
  scalnames[2] = "POrder";
  scalnames[3] = "Error";
//   scalnames[4] = "Laplac";
  TPZVec<std::string> vecnames(1);
  vecnames[0] = "Derivate";
  std::string plotfile = "singular.dx";
  an.DefineGraphMesh(2, scalnames, vecnames, plotfile);
  an.PostProcess(0);  

  delete cmesh;
  delete gmesh;
  return 0;
}

