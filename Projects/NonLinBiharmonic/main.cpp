/**
 * @param type: 1 = quadrados; 2 = triangulos
 * @param resolution: número de refinamentos uniformes
 * @param porder: ordem de interpolacao
 * @since  Fev 01, 2005
 * @author Paulo Rafael Bösing
 */
#include <pzgmesh.h>
#include <pzcmesh.h>
#include <pzcompel.h>
#include <pzgeoel.h>
#include <pzvec.h>
#include <pzgnode.h> 
#include <pzgeoelbc.h>
#include <pzpoisson3d.h>
#include <pzmatrix.h>
#include <TPZCompElDisc.h>
#include <pzfstrmatrix.h>
#include <pzstepsolver.h>
#include <fstream>
#include <pzbstrmatrix.h>
#include <pzblockdiag.h>
#include <pzbdstrmatrix.h>
#include <pzbndmat.h>
#include <pznonlinbiharmonic.h>
#include "TPZGeoElement.h"
#include "TPZInterfaceEl.h"
#include <pzanalysiserror.h>
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
#include "TPZShapeDisc.h"
#include "pzvec.h"
#include "pzquad.h"
#include "pzbndcond.h"
#include "TPZParFrontStructMatrix.h"
#include "TPZSpStructMatrix.h" 
#include "pznonlinanalysis.h" 
#include "checkconv.h"


static REAL Pi, lbd, Re;

TPZCompMesh * CreateMesh(int , int , int );

void SolveLU(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveIterative(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveIterative_2(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveIterative_GMRES(TPZNonLinearAnalysis &, TPZCompMesh *, TPZGeoMesh *);

void SolveFrontal(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha);

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] =  0.5*exp(lbd*x[0])*sin(2.*Pi*x[1])*( lbd - 2.*Pi)*(lbd + 2.*Pi)*
           (-lbd*lbd + lbd*Re + 4.*Pi*Pi)/(Re*Pi);
}

void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = x[1]-exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0); 
  deriv(0,0) = -lbd*exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0);
  deriv(1,0) = (1.0-exp(lbd*x[0])*cos(2.0* Pi*x[1]));
  deriv(2,0) = -lbd*lbd *exp(lbd*x[0])*sin(2.0* Pi*x[1])/(Pi*2.0)
               + 2.0*exp(lbd*x[0])*sin(2.0* Pi*x[1])*Pi;
 // deriv(3,0) = //D/dx do laplaciano	       
 // deriv(4,0) = //D/dy do laplaciano
  //deriv(5,0) = -lbd*exp(lbd*x[0])*cos(2.0* Pi*x[1]);// dxdy
  
}

void CC1(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 0.0;  // Cond. de Cont. de Dirichlet
  f[1] = -1. + exp(lbd*x[0]);  // Cond. de Cont. de Neumann
}


void CC2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] =  x[1]-exp(lbd*1.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
  f[1] = -lbd*exp(lbd*1.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
}
void CC3(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 2.0;  // Cond. de Cont. de Dirichlet
  f[1] = 1. - exp(lbd*x[0]);  // Cond. de Cont. de Neumann
}


void CC4(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = x[1]-exp(-lbd*0.5)*sin(2.0* Pi*x[1])/(Pi*2.0);
  f[1] = lbd*exp(-lbd*0.5)*sin(2.0*Pi*x[1])/(Pi*2.0);
}

int main(){

  TPZInterfaceElement::SetCalcStiffPenalty();

  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int ordem_interp = 3;

  cout <<" Valor de p: " << endl;
  cin >> ordem_interp;
/*  
  int numshape = (ordem_interp+1)*(ordem_interp+1);

  TPZFNMatrix<1000> phi(numshape,1),dphi(3,numshape),phix(ordem_interp+1,1),
  dphix(3,ordem_interp+1);
  TPZIntQuad oned(30,30);
  TPZManVector<REAL> X(2,0.),X0(2,0.);
  REAL weight;
  oned.Point(0,X,weight);
  REAL C = 1./16.;
  X[0] *= C;
  X[1] *= C;
  int degree = ordem_interp;
  TPZShapeDisc::MShapeType type = TPZShapeDisc::ETensorial;
  TPZShapeDisc::Shape2DFull(C,X0,X,degree,phi,dphi, type);
  TPZShapeDisc::fOrthogonal(C,X0[0],X[0],ordem_interp,phix,dphix,3);
  
  phi.Print("Funcoes de forma");
  dphi.Print("Derivadas");
  
  phix.Print("Funcao 1d");
  dphix.Print("Deriv Funcao 1d");
  return 0;
 */
  Pi = 4.0*atan(1.0);
  Re = TPZNonLinBiharmonic::Re ;
  cout<<" Reynolds = "<<Re<<endl;  
  lbd = Re/2. - sqrt(Re*Re/4. + 4.*Pi*Pi); 
  cout<<" lbd = "<<lbd<<endl;
  
  TPZGeoMesh *geomalha = NULL;
  TPZCompMesh *malha = NULL;
  malha = CreateMesh(1, 0, ordem_interp);
  geomalha = malha->Reference();

  ofstream nonlin("AnalysisNL.out");
    
  TPZNonLinearAnalysis an(malha, nonlin);
  /*TPZBandStructMatrix mat(malha);
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);
  cout << "BandStructMatrix " << endl;
  TPZStepSolver solv;
  solv.SetDirect(ELU);
  cout << "ELU " << endl;
  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  
  an.Run();
  */
  
  /*
  TPZFMatrix state(an.Solution()),range;
  range.Resize(state.Rows(),state.Cols());
  int i;
  for(i=0; i<range.Rows(); i++){
     range(i,0) = 0.01;
     state(i,0) = 0.1;
  }
  
  TPZVec<REAL> coefs(1,1.);
  
  CheckConvergence<TPZNonLinearAnalysis>(an,state,range,coefs);
*/
  //return 0;
  //SolveIterative(an, malha, geomalha);
  //SolveIterative_2(an, malha, geomalha);
  SolveIterative_GMRES(an, malha, geomalha);
  //SolveLU(an, malha, geomalha);
  an.PostProcess(1);  /* O parametro é o número de pontos p/ o "dx" construir o
			 gráfico. Ex. 1 => somente os nós, 2 => nós e pto médio */

  ofstream analy("Analysis.out");
  an.Print("Analysis.out", analy);

  ofstream geofile("GeoMalha.out");
  geomalha->Print(geofile);

  ofstream cfile("CompMalha.out");
  malha->Print(cfile);

  /* Saida de Dados */
  /* 1. Coeficientes da Soluçao Numérica */
  ofstream file("Solution.out");
  TPZFMatrix toprint = an.Solution();
  toprint.Print("solution", file);

  /* 2. Saida para dx.
        "Solution" e "Derivate" estão definidos em
	TPZMatPoisson3d::VariableIndex  */
  TPZVec<char *> scalar_names(1);
  TPZVec<char *> vec_names(1);
  scalar_names[0] = "Solution";
  vec_names[0] = "Derivate";
  an.DefineGraphMesh(2, scalar_names, vec_names, "Solution.dx");

  an.PostProcess(2);  /* O parametro é o número de pontos p/ o "dx" construir o
			 gráfico. Ex. 1 => somente os nós, 2 => nós e pto médio */

  ofstream out("teste.dat");


  /* Plotar as normas dos erros */
  an.SetExact(ExactSolution);
  TPZVec<REAL> pos;
  an.PostProcess(pos,out);

  delete malha;
  delete geomalha;

}//end of main


TPZCompMesh * CreateMesh(int type, int resolution, int porder){

  TPZCompElDisc::gDegree = porder;

  TPZGeoMesh *geomalha = new TPZGeoMesh;

  TPZVec<REAL> coord(2);

 int nodind;

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = -0.5;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(0, coord, *geomalha);

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.5;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(1, coord, *geomalha);


  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = -0.5;
  coord[1] = 2.0;
  geomalha->NodeVec()[nodind].Initialize(2, coord, *geomalha);


  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.5;
  coord[1] = 2.0;
  geomalha->NodeVec()[nodind].Initialize(3, coord, *geomalha);

  TPZVec<TPZGeoEl *> geoel(1);
  int index;
  TPZVec<int> incid(4);
  incid[0] = 0;
  incid[1] = 1;
  incid[2] = 3;
  incid[3] = 2;
  geoel[0] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);


  geomalha->BuildConnectivity();

  //   Refinar malha
    TPZVec<TPZGeoEl *> children, netos, bisnetos,bisbisnetos,bisbisbisnetos;


     for(int i=0;i<1;i++){
       geoel[i]->Divide(children);
         for(int j = 0; j < children.NElements(); j++){
            children[j]->Divide(netos);
                  for(int j2 = 0; j2< netos.NElements(); j2++) {/*
                	 netos[j2]->Divide(bisnetos);
                          for(int j3 = 0; j3< bisnetos.NElements(); j3++) {
               	           bisnetos[j3]->Divide(bisbisnetos);
    			     for(int j4 = 0; j4< bisbisnetos.NElements(); j4++) {
               	               bisbisnetos[j4]->Divide(bisbisbisnetos);
    			     }
         		 }*/
		 }
    
    }
   }
 

  TPZGeoElBC elbc1(geoel[0],4,-1,*geomalha);
  TPZGeoElBC elbc2(geoel[0],5,-2,*geomalha);
 
  TPZGeoElBC elbc3(geoel[0],6,-3,*geomalha);
  TPZGeoElBC elbc4(geoel[0],7,-4,*geomalha);

//   TPZGeoElBC elbc5(geoel[2],5,-2,*geomalha);
//   TPZGeoElBC elbc6(geoel[2],6,-3,*geomalha);

//   TPZGeoElBC elbc7(geoel[3],6,-3,*geomalha);
//   TPZGeoElBC elbc8(geoel[3],7,-2,*geomalha);

  TPZCompMesh * malha = new TPZCompMesh(geomalha);
  malha->SetDimModel(2); 

  TPZNonLinBiharmonic *mater;
  mater = new TPZNonLinBiharmonic(1,0.);  // segundo par. é a f(x)
                                   // primeiro par. é o material

  mater->SetForcingFunction(Forcing1);
  malha->InsertMaterialObject(mater);

 TPZFMatrix val1(1,1,0.), val2(2,1,0.);
  
 TPZMaterial *bnd1 = mater->CreateBC(-1,0, val1, val2);
 TPZMaterial *bnd2 = mater->CreateBC(-2,0, val1, val2);
 TPZMaterial *bnd3 = mater->CreateBC(-3,0, val1, val2);
 TPZMaterial *bnd4 = mater->CreateBC(-4,0, val1, val2);
 
 bnd1->SetForcingFunction(CC1);
 malha->InsertMaterialObject(bnd1);
 bnd2->SetForcingFunction(CC2);
 malha->InsertMaterialObject(bnd2);
 bnd3->SetForcingFunction(CC3);
 malha->InsertMaterialObject(bnd3);
 bnd4->SetForcingFunction(CC4);
 malha->InsertMaterialObject(bnd4);

/*  
  TPZBndCond *cond_front[4];

  cond_front[0] = mater->CreateBC(-1, 0, val1, val2);
  cond_front[1] = mater->CreateBC(-2, 0, val1, val2);
  cond_front[2] = mater->CreateBC(-3, 0, val1, val2);
  cond_front[3] = mater->CreateBC(-4, 0, val1, val2);

  cond_front[0]->SetForcingFunction(CC1);
  cond_front[1]->SetForcingFunction(CC2); 
  cond_front[2]->SetForcingFunction(CC3);
  cond_front[3]->SetForcingFunction(CC4);


//                                             // 1 cond. fronteira (negativo)
//                                             // 2 Dirichlet
//                                             // 3 Caso as condicoes sáo mistas
//                                             // 4 Dirichlet ou Newmann
//   // Se a cond. de front. for uma funcao , ver TPZMaterial::fForcingFunction para BndCond::Contribute
  // cond_front[0]->SetForcingFunction(Dirichlet1);             Tem algun problema.
  //  cond_front[1]->SetForcingFunction(Dirichlet2);

  malha->InsertMaterialObject(mater);
  // malha->InsertMaterialObject(cond_front[0]);
  malha->InsertMaterialObject(cond_front[0]);
  malha->InsertMaterialObject(cond_front[1]);
  malha->InsertMaterialObject(cond_front[2]);
  malha->InsertMaterialObject(cond_front[3]);
*/
  
  /* Para usar elementos descontinuos. Caso for subtraido, será elementos continuos. */
   TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
   TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);


  malha->AutoBuild();
  
  return malha;

}//end of CreateMesh

void SolveLU(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){
  TPZBandStructMatrix mat(malha);
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);
  cout << "BandStructMatrix " << endl;
  TPZStepSolver solv;
  solv.SetDirect(ELU);
  cout << "ELU " << endl;
  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  cout << endl;
  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;
  
  ofstream iterproc("IteraProc.dat");
  
  
  an.IterativeProcess(iterproc,0.0000001,50);
  //an.Assemble();  /* an.Run() faz a Assemble e Solve */
  //an.Solve();
  cout << endl;
  /*TPZStepSolver solv2;
  TPZSpStructMatrix sparse(malha);
  an.SetStructuralMatrix(sparse);
  solv2.SetGMRES( 2000, 10, solv, 1.e-16, 0);
  an.SetSolver(solv2);
  an.Assemble();
  cout << "Solve " << endl;
  an.Solve();
*/
  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  TPZBandStructMatrix mat(malha)";
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);

    //     bandmat->Print();
 
  }

}

void SolveIterative(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){

  TPZBandStructMatrix mat(malha);
  //  TPZFStructMatrix mat(malha);
  an.SetStructuralMatrix(mat);
  cout << "BandStructMatrix " << endl;

  TPZBlockDiagonalStructMatrix strblock(malha);
  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strblock.AssembleBlockDiagonal(*block);

  TPZStepSolver solv, precond(block);

  precond.SetDirect(ELU);  
  //TPZBandStructMatrix mat(malha);
  //TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);
  //TPZParFrontStructMatrix mat(malha);
  cout << "BandStructMatrix " << endl;
  
  cout << " GMRES " << endl;

  solv.SetGMRES( 2000, 10, precond, 0.0000000001, 0);

  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  cout << endl;
  cout << "Assemble " << endl;

  //an.Assemble();  /* an.Run() faz a Assemble e Solve */
  cout << endl;
  cout << "Solve " << endl;

  //an.Solve();

  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  = " << bandmat->Rows() << endl;

    //     bandmat->Print();
  }

}

void SolveIterative_2(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){

  TPZParFrontStructMatrix<TPZFrontNonSym> mat(malha);

  cout << " TPZParFrontStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  cout << endl;
  
  TPZMatrix *sec_mat = NULL;
  
  TPZStepSolver solv;
  solv.SetDirect(ELU);
  cout << "ELU " << endl;
  solv.SetMatrix(sec_mat);
  an.SetSolver(solv);
  
  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;
  
  ofstream iterproc("IteraProc.dat");
  
  
  an.IterativeProcess(iterproc,0.0000001,50);
  cout << endl;
 
  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  TPZBandStructMatrix mat(malha)";
  }
  
}

void SolveIterative_GMRES(TPZNonLinearAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){

 
  TPZSpStructMatrix mat(malha);
 
  cout << " TPZSpStructMatrix " << endl;
  an.SetStructuralMatrix(mat);
  cout << endl;
  
  TPZMatrix *sec_mat = NULL;
  
  TPZStepSolver Pre;
  
  TPZBlockDiagonalStructMatrix strBlockDiag(malha);
  
  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strBlockDiag.AssembleBlockDiagonal(*block);
  
  Pre.SetMatrix(block);
  Pre.SetDirect(ELU);
  
  TPZStepSolver solv;
  solv.SetGMRES(3000,200,Pre,1e-9,0);
  cout << "GMRES " << endl;
  solv.SetMatrix(sec_mat);
  an.SetSolver(solv);
  
  an.Solution().Redim(0,0);
  cout << "Assemble " << endl;
  
  ofstream iterproc("IteraProc.dat");
  for(int j1=1; j1<6; j1++){
    TPZNonLinBiharmonic::Re = j1*200;
    an.IterativeProcess(iterproc,0.0000001,50);
  } 
  cout << endl;
 
  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  TPZBandStructMatrix mat(malha)";
  }

  
}
 