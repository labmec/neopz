//$Id: main.cpp,v 1.2 2005-02-04 19:47:07 tiago Exp $

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

static REAL epsilon;

void ForcingFunction(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  REAL x = pto[0];
  REAL y = pto[1];
  
  force.Resize(1);
  force[0] = (-(-1.0 + exp(1.0/epsilon)) * epsilon * (-2.0+x+y)
             +exp( (x+y-x*y)/epsilon ) * (-x + x*x + (-1+y)*y) ) / ( ( -1.0 + exp(1.0/epsilon) ) * epsilon);
  force[0] *= -1.0;
}

void ForcingFunctionSemBeta(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  REAL x = pto[0];
  REAL y = pto[1];
  
  force.Resize(1);
  force[0] = exp((x+y-x*y)/epsilon)*(2.-2.*x+x*x-2.*y+y*y)/( (-1.+exp(1./epsilon))*epsilon );
}

void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  deriv.Resize(2,1);
  
  u[0] = ( exp(-1.0/epsilon) -exp( (-1.0+x) * (1.0-y) / epsilon ) ) / ( 1.0 - exp( -1.0/epsilon) ) + x + (1.0-x)*y;
  
  deriv[0] = ( exp((x+y-x*y)/epsilon) + epsilon - exp(1.0/epsilon) * epsilon) * (-1.0 + y) / ( (-1.0+exp(1.0/epsilon))*epsilon );
  
  deriv[1] = ( exp((x+y-x*y)/epsilon) + epsilon - exp(1.0/epsilon) * epsilon) * (-1.0 + x) / ( (-1.0+exp(1.0/epsilon))*epsilon );
  
}

void Dirichlet_X_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u) {
//  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  u[0] = -(-1.0+exp(y/epsilon)+y-exp(1.0/epsilon)*y) / (-1.0 + exp(1.0/epsilon));
}

void Dirichlet_Y_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u) {
  REAL x = pto[0];
//  REAL y = pto[1];
  u.Resize(1);
  u[0] = -(-1.0+exp(x/epsilon)+x-exp(1.0/epsilon)*x) / (-1.0 + exp(1.0/epsilon));
}

//4 quadrados iguais, podendo ser divididos
TPZCompMesh *CreateMesh(int h = 0);

//8 triangulos
TPZCompMesh *CreateMesh2();

////////////////////////////////////////////////////////////////////////////////


int main(){
   
   int nmaxeq = 10000;//numero maximo de equacoes do problema
for(int pp = 0; pp < 8; pp++){
 for(int hh = 0; hh < 6; hh++){

  TPZCompElDisc::SetOrthogonalFunction(TPZShapeDisc::Legendre);

  int p, h;
  p = pp;
  epsilon = 0.01;
  h = hh;
      
  TPZCompEl::gOrder = p;
  TPZCompElDisc::gDegree = p;

  TPZCompMesh *cmesh;
  cmesh = CreateMesh(h); //    cmesh = CreateMesh2();
    
  TPZGeoMesh *gmesh = cmesh->Reference();

  TPZAnalysis an(cmesh);

  TPZParFrontStructMatrix <TPZFrontNonSym> full(cmesh);
//  full.SetNumberOfThreads(3);
//  TPZFStructMatrix full(cmesh);  
#define direct
#ifdef direct
//  TPZFStructMatrix/*TPZSkylineStructMatrix*/ full(cmesh);
  an.SetStructuralMatrix(full);
  TPZStepSolver step;
  step.SetDirect( ELU /*ECholesky*/ /*ELDLt*/ );
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
  step.SetGMRES( 200000, 20, precond, 1.e-13, 0 ); 
  //step.SetCG( 200000, *precond, 1e-13, 0);
  //delete precond;
  an.SetSolver(step);
#endif

  cout << "\nNumero de equacoes: " << an.Mesh()->NEquations() << endl;
   if (an.Mesh()->NEquations() > nmaxeq) { 
     cout << "skipping simulation...\n" << endl;
     break;
   }
  char filename[20];
  sprintf(filename,"baumann_p%d_h%d.dat",p,h);
  char filedx[20];
  sprintf(filedx,"baumann_p%d_h%d.dx",p,h);

  
  
/*{  
  TPZFMatrix fillin;
  cmesh->ComputeFillIn(50,fillin);
  //fillin.Print("Fillin of the computable mesh");
  VisualMatrix(fillin , filedx);
}*/

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

TPZCompMesh *CreateMesh(int h) {

  REAL co[9][2] = {{0.5,0.5},{0.5,0.},{1.,0.},{1.,0.5},{1.,1.},{0.5,1.},{0.,1.},{0.,0.5},{0.,0.}};
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

  TPZVec<TPZGeoEl*> filhos;
  for(int i = 0; i < h; i++){
     int n = gmesh->NElements();
     for(int j = 0; j < n; j++){
	if (gmesh->ElementVec()[j]->Dimension() == 2) gmesh->ElementVec()[j]->Divide(filhos);
     }
  }

  //TPZGeoElBC gbc;

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-1,*gmesh); // bottom
    
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMatPoisson3d *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(ForcingFunction);
  mat->SetNonSymmetric();

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;

  REAL beta = 1.0;
  mat->SetParameters(epsilon, beta, convdir);
  int nstate = 1;
  
  
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[3];
  val2.Zero();
  bc[0] = mat->CreateBC(-3, 0,val1,val2);
  
  bc[1] = mat->CreateBC(-1, 0, val1, val2);
  bc[1]->SetForcingFunction(Dirichlet_Y_IgualA_0);
  
  bc[2] = mat->CreateBC(-2, 0, val1, val2);
  bc[2]->SetForcingFunction(Dirichlet_X_IgualA_0);
  
  cmesh->InsertMaterialObject(mat);
  for(int ii = 0; ii < 3; ii++) cmesh->InsertMaterialObject(bc[ii]);

  
     TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
     TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc); 

  
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
  mat->SetForcingFunction(ForcingFunction);

  TPZManVector<REAL,2> convdir(2,0.);
  convdir[0] = 1.;
  convdir[1] = 1.;
  mat->SetParameters(epsilon, 0., convdir);
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


