//$Id: main.cpp,v 1.1 2005-02-03 18:07:23 tiago Exp $

/**
 * Galerkin descontinuo: problema de camada limite
 * 03/02/2005
 */


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

#include <time.h>
#include <stdio.h>


void ForcingFunction(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = - exp(0.75 * (x[0] + x[1])) * (8. * (1. - x[1] * x[1]) + 12. * x[0] * (1. - x[1] * x[1]) -4.5 * (1. - x[0] * x[0])* (1. - x[1] * x[1] )+
					  8. * (1. - x[0] * x[0]) + 12. * x[1] * (1. - x[0] * x[0]));
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = 4. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) * exp(0.75 * (x[0] + x[1]));
  deriv(0,0) = (-8. * x[0] * (1. - x[1] * x[1]) + 3. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) ) * exp(0.75 * (x[0] + x[1]));
  deriv(1,0) = (-8. * x[1] * (1. - x[0] * x[0]) + 3. * (1. - x[0] * x[0] ) * (1. - x[1] * x[1]) ) * exp(0.75 * (x[0] + x[1]));
}

void Dirichlet1(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = exp(-(x[0]*x[0]+1));
}

void Dirichlet2(TPZVec<REAL> &x, TPZVec<REAL> &f) {
  f[0] = exp(-(x[1]*x[1]+1));
}

//4 quadrados iguais, podendo ser divididos
TPZCompMesh *CreateMesh();

//8 triangulos
TPZCompMesh *CreateMesh2();

static REAL PI;

////////////////////////////////////////////////////////////////////////////////


int main(){
   
  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int p, h;
  char filename[20];
  char filedx[20];
  PI = 4.*atan(1.);

  ofstream out(filename);

  TPZCompEl::gOrder = p;
  TPZCompElDisc::gDegree = p;

  TPZCompMesh *cmesh;
  cmesh = CreateMesh(); //    cmesh = CreateMesh2();
    
  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZAnalysis an(cmesh);

//  TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
//  full.SetNumberOfThreads(3);
//  TPZFStructMatrix full(cmesh);  
#define direct
#ifdef direct
  TPZBandStructMatrix/*TPZSkylineStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect( ELU /*ELDLt*/ );
  an.SetSolver(step);
#endif

//#define iter
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  TPZSpStructMatrix full(cmesh);
  //TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);
  
//  TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , true);
  TPZCopySolve precond( full.Create() );
  step.ShareMatrix( precond );  
  step.SetGMRES( 200000, 100, precond, 1.e-13, 0 ); 
  //step.SetCG( 200000, *precond, 1e-13, 0);
  //delete precond;
  an.SetSolver(step);
#endif

  cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
   if (an.Mesh()->NEquations() > 22000) { 
     cout << "skipping simulation...\n" << endl;
     continue;
   }
  char filename[20];
  sprintf(filename,"baumann_p%d_h%d.dat",p,h);
  char filedx[20];
  sprintf(filedx,"baumann_p%d_h%d.dx",p,h);
  ofstream out(filename);
  
  
{  
  TPZFMatrix fillin;
  cmesh->ComputeFillIn(50,fillin);
  //fillin.Print("Fillin of the computable mesh");
  VisualMatrix(fillin , filedx);
}

  an.Run();

/**** Aqui faz DX ****/
/*  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);

  if (gMeshType == 1)
    an.PostProcess(2);

  if (gMeshType == 2)
    an.PostProcess(0);
*/

  an.SetExact(ExactSolution);
  TPZVec<REAL> pos;
  an.PostProcess(pos,out);
  out << "\nNumero de equacoes: " << an.Solution().Rows() << endl;  
  
  delete cmesh;
  delete gmesh;
  return 0;
}




//////////////////////// COMECAM AS MALHAS:
/////////////////////// CREATEMESH -> 4 QUADRADOS
////////////////////// CREATEMESH2 -> 8 TRIANGULOS



#include "TPZGeoElement.h"

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

TPZCompMesh *CreateMesh() {

  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int indices[4][4] = {{1,2,3,0},{0,3,4,5},{7,0,5,6},{8,1,0,7}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = 4;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<nelem; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod]=indices[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index);
  }

  gmesh->BuildConnectivity();

  TPZGeoElBC gbc;

  TPZGeoElBC gbc1(elvec[0],4,-3,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-3,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-3,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-3,*gmesh); // bottom
  

  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);
  mat->SetNonSymmetric();


  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  //Difusao pura. conveccao = 0.
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc;
  
  val2.Zero();
  bc = mat->CreateBC(-3, 0,val1,val2);

  
  cmesh->InsertMaterialObject(mat);
  cmesh->InsertMaterialObject(bc);



  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc); 
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  cout << "\nGALERKIN DESCONTINUO" << endl;
//  cout << "\nELEMENTOS FINITOS" << endl;
  
  cmesh->AutoBuild();
  //  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}


TPZCompMesh *CreateMesh2() {
  REAL co[9][2] = {{0.,0.},{0.,-1.},{1.,-1.},{1.,0.},{1.,1.},{0.,1.},{-1.,1.},{-1.,0.},{-1,-1}};
  int numquad = 0;
  int indicesquad[2][4] = {{0,1,2,3},{0,5,6,7}};
  int numtriang = 8;
  int indicestriang[8][3] = {{0,2,3},{0,1,2},{0,3,5},{3,4,5},{0,7,1},{7,8,1},{0,5,6},{0,6,7}};
  TPZGeoEl *elvec[4];
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  int nnode = 9;
  int nelem = numquad + numtriang;
  int nod;
  for(nod=0; nod<nnode; nod++) {
    int nodind = gmesh->NodeVec().AllocateNewElement();
    TPZManVector<REAL,2> coord(2);
    coord[0] = co[nod][0];
    coord[1] = co[nod][1];
    gmesh->NodeVec()[nodind].Initialize(nod,coord,*gmesh);
  }

  int el;
  for(el=0; el<numquad; el++) {
    TPZManVector<int,4> nodind(4);
    for(nod=0; nod<4; nod++) nodind[nod] = indicesquad[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(EQuadrilateral,nodind,1,index); 
  }
 
  for(el=0; el<numtriang; el++) {
    TPZManVector<int,3> nodind(3);
    for(nod=0; nod<3; nod++) nodind[nod]=indicestriang[el][nod];
    int index;
    elvec[el] = gmesh->CreateGeoElement(ETriangle,nodind,1,index);
  }


  gmesh->BuildConnectivity();




  TPZGeoElBC gbc1(elvec[1],4,-3,*gmesh); 
  TPZGeoElBC gbc2(elvec[0],4,-3,*gmesh); 
  TPZGeoElBC gbc3(elvec[3],3,-3,*gmesh); 
  TPZGeoElBC gbc4(elvec[3],4,-3,*gmesh); 
  TPZGeoElBC gbc5(elvec[6],4,-3,*gmesh); 
  TPZGeoElBC gbc6(elvec[7],4,-3,*gmesh); 
  TPZGeoElBC gbc7(elvec[5],3,-3,*gmesh); 
  TPZGeoElBC gbc8(elvec[5],4,-3,*gmesh); 

  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(gDif, 0., convdir);
  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);

  TPZBndCond *bc;
  val2.Zero(); val1.Zero(); bc = mat->CreateBC(-3, 0, val1, val2);
  
  cmesh->InsertMaterialObject(mat);

  cmesh->InsertMaterialObject(bc);

  
  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  
  return cmesh;
}


int main(){

TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);
  
cout << "Script" << endl;

for(int tiagop = 2; tiagop < 9; tiagop++){
   
  for( int tiagoh = 0; tiagoh < 7; tiagoh++){

  int p, h;
  p = tiagop;
  h = tiagoh;
  
  cout << "\n\n" << "p = " << p << " - h = " << h << endl;
  
 //cout << "\nQuadrado = 1; Triang = 2" << endl;
  gMeshType = 1;

 gDif = 1.0;

  if (h == 0) {
    gDivide[0] = 0; 
    gDivide[1] = 0;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 1) {
    gDivide[0] = 1;
    gDivide[1] = 0;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  } 
  if (h == 2) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 0;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }  
  if (h == 3) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 0;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 4) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 0;
    gDivide[5] = 0;
  }
  if (h == 5) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 1;
    gDivide[5] = 0;
  }
  if (h == 6) {
    gDivide[0] = 1;
    gDivide[1] = 1;
    gDivide[2] = 1;
    gDivide[3] = 1;
    gDivide[4] = 1;
    gDivide[5] = 1;
  }


  PI = 4.*atan(1.);

  TPZCompEl::gOrder = p;
  TPZCompElDisc::gDegree = p;

  TPZCompMesh *cmesh;
  cmesh = CreateMesh();

  TPZGeoMesh *gmesh = cmesh->Reference();
  
//  gmesh->Print();

  TPZAnalysis an(cmesh);


//  TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
//  full.SetNumberOfThreads(3);
//  TPZFStructMatrix full(cmesh);  
#define direct
#ifdef direct
  TPZBandStructMatrix/*TPZSkylineStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect( ELU /*ELDLt*/ );
  an.SetSolver(step);
#endif

//#define iter
#ifdef iter
  cout << "ITER_SOLVER" << endl;  
  TPZSpStructMatrix full(cmesh);
  //TPZSkylineStructMatrix full(cmesh);
  an.SetStructuralMatrix(full);  
  TPZStepSolver step( full.Create() );
  an.SetSolver(step);
  
//  TPZMatrixSolver * precond = an.BuildPreconditioner(TPZAnalysis::EElement , true);
  TPZCopySolve precond( full.Create() );
  step.ShareMatrix( precond );  
  step.SetGMRES( 200000, 100, precond, 1.e-13, 0 ); 
  //step.SetCG( 200000, *precond, 1e-13, 0);
  //delete precond;
  an.SetSolver(step);
#endif

  cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
   if (an.Mesh()->NEquations() > 22000) { 
     cout << "skipping simulation...\n" << endl;
     continue;
   }
  char filename[20];
  sprintf(filename,"baumann_p%d_h%d.dat",p,h);
  char filedx[20];
  sprintf(filedx,"baumann_p%d_h%d.dx",p,h);
  ofstream out(filename);

{  
  TPZFMatrix fillin;
  cmesh->ComputeFillIn(50,fillin);
  //fillin.Print("Fillin of the computable mesh");
  VisualMatrix(fillin , filedx);
}

 // PIRAArrumaNormal( *cmesh );
 
  an.Run();

  an.SetExact(ExactSolution);
  TPZVec<REAL> pos;
  an.PostProcess(pos,out);
  out << "\nNumero de equacoes: " << an.Solution().Rows() << endl;  
  
  delete cmesh;
  delete gmesh;
  }//p
  }//h
  return 0;
}

#include "TPZInterfaceEl.h"
void PIRAArrumaNormal( TPZCompMesh &cmesh ){

  cout << "PIRAArrumaNormal" << endl;
  int nel = cmesh.ElementVec().NElements();
  for(int i = 0; i < nel; i++){
    TPZInterfaceElement * face = dynamic_cast<TPZInterfaceElement*> (cmesh.ElementVec()[i]);
    if (face) face->PIRAArrumaNormal();    
  }







}