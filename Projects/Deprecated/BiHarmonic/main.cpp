/**
 * @param type: 1 = quadrados; 2 = triangulos
 * @param resolution: número de refinamentos uniformes
 * @param porder: ordem de interpolacao
 * @since Mar 22, 2004
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
#include <pzbiharmonic.h>
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


static REAL Pi;

TPZCompMesh * CreateMesh(int type, int resolution, int porder);

void SolveLU(TPZAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha);

void SolveIterative(TPZAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha);

void SolveFrontal(TPZAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha);

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = -4.*Pi*Pi*Pi*Pi*( cos(2.*Pi*x[0]) - 2.*cos( 2.*Pi*(x[0] - x[1]) ) +
	                      cos(2.*Pi*x[1]) - 2.*cos( 2.*Pi*(x[0] + x[1]) ) );
}

void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = sin(Pi*x[0])*sin(Pi*x[0])*sin(Pi*x[1])*sin(Pi*x[1]);
  //  deriv.Redim(3,1);
  deriv(0,0) = 2.*Pi*sin(Pi*x[0])*cos(Pi*x[0])*sin(Pi*x[1])*sin(Pi*x[1]);
  deriv(1,0) = 2.*Pi*sin(Pi*x[1])*cos(Pi*x[1])*sin(Pi*x[0])*sin(Pi*x[0]);
  deriv(2,0) = 2.*Pi*Pi*cos(Pi*x[1])*cos(Pi*x[1])*sin(Pi*x[0])*sin(Pi*x[0]) +
               2.*Pi*Pi*cos(Pi*x[0])*cos(Pi*x[0])*sin(Pi*x[1])*sin(Pi*x[1]) -
               4.*Pi*Pi*sin(Pi*x[0])*sin(Pi*x[0])*sin(Pi*x[1])*sin(Pi*x[1]);	  
}

void CC1(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 0.;  // Cond. de Cont. de Dirichlet
  f[1] = 0.;  // Cond. de Cont. de Neumann
}


void CC2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = 0.; //x[0]*(-1. + x[0]);
  f[1] = 0.;
}

int main2(){

  Pi = 4.*atan(1.);
 
  TPZVec<REAL> ponto(2);
  TPZVec<REAL> u(1);
  TPZFMatrix deriv(3,1);
  for (int i = 0; i < 11; i++)
    for (int j = 0; j < 11; j++){
      ponto[0] = (double) 0.1 * i;
      ponto[1] = (double) 0.1 * j;
      ExactSolution(ponto, u, deriv);
      cout << "\nPonto = " << ponto[0] << " , " << ponto[1] << endl;
      cout << "u = " << u[0] << endl;
      cout << "deriv = " << deriv(0,0) <<   " , "  << deriv(1,0) <<  " , " << deriv(2,0) << endl;
    }
}

int main(){
  
   TPZInterfaceElement::SetCalcStiffPenalty();

   TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int ordem_interp = 3;

  cin >> ordem_interp;

  Pi = 4.*atan(1.);
 
  TPZGeoMesh *geomalha = NULL; 
  TPZCompMesh *malha = NULL;
  malha = CreateMesh(1, 0, ordem_interp);
  geomalha = malha->Reference();
 
  TPZAnalysis an(malha);

  //  SolveIterative(an, malha, geomalha);
  SolveLU(an, malha, geomalha);
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
  an.DefineGraphMesh(2, scalar_names, vec_names, "Solutionp6ef.dx");
  
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
  coord[0] = 0.;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(0, coord, *geomalha);
  
//    nodind  =  geomalha->NodeVec().AllocateNewElement();
//    coord[0] = 0.5;
//    coord[1] = 0.;
//    geomalha->NodeVec()[nodind].Initialize(1, coord, *geomalha);

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.;
  coord[1] = 0.;
  geomalha->NodeVec()[nodind].Initialize(1, coord, *geomalha); 

//   nodind  =  geomalha->NodeVec().AllocateNewElement();
//   coord[0] = 0.;
//   coord[1] = 0.5;
//   geomalha->NodeVec()[nodind].Initialize(3, coord, *geomalha);

//   nodind  =  geomalha->NodeVec().AllocateNewElement();
//   coord[0] = 0.5;
//   coord[1] = 0.5;
//   geomalha->NodeVec()[nodind].Initialize(4, coord, *geomalha);

//   nodind  =  geomalha->NodeVec().AllocateNewElement();
//   coord[0] = 1.;
//   coord[1] = 0.5;
//   geomalha->NodeVec()[nodind].Initialize(5, coord, *geomalha);

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 0.;
  coord[1] = 1.;
  geomalha->NodeVec()[nodind].Initialize(2, coord, *geomalha);

//    nodind  =  geomalha->NodeVec().AllocateNewElement();
//    coord[0] = 0.5;
//    coord[1] = 1.;
//    geomalha->NodeVec()[nodind].Initialize(7, coord, *geomalha);

  nodind  =  geomalha->NodeVec().AllocateNewElement();
  coord[0] = 1.;
  coord[1] = 1.;
  geomalha->NodeVec()[nodind].Initialize(3, coord, *geomalha);

  TPZVec<TPZGeoEl *> geoel(1);
  int index;
  TPZVec<int> incid(4);
  incid[0] = 0;
  incid[1] = 1;
  incid[2] = 3;
  incid[3] = 2;
  geoel[0] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);

//   incid[0] = 1;
//   incid[1] = 2;
//   incid[2] = 5;
//   incid[3] = 4;
//   geoel[1] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);

//   incid[0] = 4;
//   incid[1] = 5;
//   incid[2] = 8;
//   incid[3] = 7;
//   geoel[2] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);

//   incid[0] = 3;
//   incid[1] = 4;
//   incid[2] = 7;
//   incid[3] = 6;
//   geoel[3] = geomalha->CreateGeoElement(EQuadrilateral, incid, 1, index);


  geomalha->BuildConnectivity();

//   TPZVec<TPZGeoEl*> children;
//   geoel[0]->Divide(children);

//   for (int ii = 0; ii < children.NElements(); ii++){
//     for (int inode = 0; inode < children[ii]->NNodes(); inode++)
//       if (children[ii]->NodePtr(inode)->Id() == 0){
// 	TPZVec<TPZGeoEl*> children2;
// 	children[ii]->Divide(children2);
// 	for(int ich = 0; ich < children2.NElements(); ich++)
// 	  for(int ichnode = 0; ichnode < children2[ich]->NNodes(); ichnode++)
// 	    if (children2[ich]->NodePtr(ichnode)->Id()==0){
// 	      TPZVec<TPZGeoEl*> children3;
// 	      children2[ich]->Divide(children3);
// 	    }
//       }//if Id==0
//   }//ii

//   geoel[2]->Divide(children);
//   for (int ii = 0; ii < children.NElements(); ii++){
//     for (int inode = 0; inode < children[ii]->NNodes(); inode++)
//       if (children[ii]->NodePtr(inode)->Id() == 8){
// 	TPZVec<TPZGeoEl*> children2;
// 	children[ii]->Divide(children2);
// 	for(int ich = 0; ich < children2.NElements(); ich++)
// 	  for(int ichnode = 0; ichnode < children2[ich]->NNodes(); ichnode++)
// 	    if (children2[ich]->NodePtr(ichnode)->Id()==8){
// 	      TPZVec<TPZGeoEl*> children3;
// 	      children2[ich]->Divide(children3);
// 	    }

//       }//if Id==0
//   }//ii


//  Refinar malha 
   TPZVec<TPZGeoEl *> children, netos;    
   geoel[0]->Divide(children);
//   for(int i = 0; i < children.NElements(); i++)
   //    children[i]->Divide(netos);

  /*  for(int i = 0; i< 4; i++) {    
    geoel[i]->Divide(children);    
    for(int j = 0; j< children.NElements(); j++) {    
      children[j]->Divide(netos);    
    }   
  } */  
  



  TPZGeoElBC elbc1(geoel[0],4,-2,*geomalha);
  TPZGeoElBC elbc2(geoel[0],5,-2,*geomalha);
 
  TPZGeoElBC elbc3(geoel[0],6,-3,*geomalha);
  TPZGeoElBC elbc4(geoel[0],7,-3,*geomalha);

//   TPZGeoElBC elbc5(geoel[2],5,-2,*geomalha);
//   TPZGeoElBC elbc6(geoel[2],6,-3,*geomalha);

//   TPZGeoElBC elbc7(geoel[3],6,-3,*geomalha);
//   TPZGeoElBC elbc8(geoel[3],7,-2,*geomalha);

  TPZCompMesh * malha = new TPZCompMesh(geomalha);
  malha->SetDimModel(2);

  TPZBiharmonic *mater;
  mater = new TPZBiharmonic(1,0.);  // segundo par. é a f(x)
                                   // primeiro par. é o material

  mater->SetForcingFunction(Forcing1);


  TPZBndCond *cond_front[3];

  TPZFMatrix val1(1,1,0.), val2(2,1,0.);

  //  cond_front[0] = mater->CreateBC(-1, 0, val1, val2);
  cond_front[1] = mater->CreateBC(-2, 0, val1, val2);
  cond_front[2] = mater->CreateBC(-3, 0, val1, val2);

  //  cond_front[0]->SetForcingFunction(Dirichlet1);             

  cond_front[1]->SetForcingFunction(CC1);
 
  cond_front[2]->SetForcingFunction(CC2);



//                                             // 1 cond. fronteira (negativo)
//                                             // 2 Dirichlet
//                                             // 3 Caso as condicoes sáo mistas
//                                             // 4 Dirichlet ou Newmann 
//   // Se a cond. de front. for uma funcao , ver TPZMaterial::fForcingFunction para BndCond::Contribute 
  // cond_front[0]->SetForcingFunction(Dirichlet1);             Tem algun problema.
  //  cond_front[1]->SetForcingFunction(Dirichlet2);
  



  malha->InsertMaterialObject(mater);
  // malha->InsertMaterialObject(cond_front[0]);
  malha->InsertMaterialObject(cond_front[1]);
  malha->InsertMaterialObject(cond_front[2]);



 
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

void SolveLU(TPZAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){
  //  TPZBandStructMatrix mat(malha);
  TPZFStructMatrix mat(malha);
  TPZStepSolver solv;
  solv.SetDirect(ELU);
  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  an.Assemble();  /* an.Run() faz a Assemble e Solve */
  an.Solve();  
}

void SolveIterative(TPZAnalysis &an, TPZCompMesh *malha, TPZGeoMesh *geomalha){
  
  TPZBandStructMatrix mat(malha);
  //  TPZFStructMatrix mat(malha);
  an.SetStructuralMatrix(mat);
  
  TPZBlockDiagonalStructMatrix strblock(malha);
  TPZBlockDiagonal * block = new TPZBlockDiagonal();
  strblock.AssembleBlockDiagonal(*block);
  
  TPZStepSolver solv, precond(block);
  
  precond.SetDirect(ELU);

  solv.SetGMRES( 2000, 10, precond, 0.0000000001, 0); 

  an.SetSolver(solv);
  an.SetStructuralMatrix(mat);
  an.Assemble();  /* an.Run() faz a Assemble e Solve */
  an.Solve();  

  TPZFBMatrix * bandmat = dynamic_cast<TPZFBMatrix*>(solv.Matrix());
  if (bandmat){
    cout << endl << "Banda = " << bandmat->GetBand() << endl;
    cout << "No equacoes = " << malha->NEquations() << endl;
    cout << "No equacoes 2  = " << bandmat->Rows() << endl;
    
    //     bandmat->Print();
  }
  
}

