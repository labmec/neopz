//$Id: main.cpp,v 1.4 2003-12-01 14:53:55 tiago Exp $
/**
 * Galerkin descontinuo: visita do professor Igor.
 * 24/11/2003
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
#include "pzbstrmatrix.h"
#include "pzstepsolver.h"    		
#include "pzonedref.h"

#include "pzpoisson3d.h"

#include <time.h>
#include <stdio.h>

static REAL PI;

int gPrintLevel = 0;

int gDivide[6];

/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] = -2.0 * PI * PI * sin( PI * x[0]) * sin(PI * x[1]) ;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = sin(PI*x[0])*sin(PI*x[1]);
  deriv(0,0) = PI*cos(PI*x[0])*sin(PI*x[1]);
  deriv(1,0) = PI*cos(PI*x[1])*sin(PI*x[0]); 
}*/

void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
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

/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  disp[0] =  PI * PI * sin(PI * (x[1]+1.)/2.) / 4.;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = sin(PI*(x[1]+1.)/2.);
  deriv(0,0) = 0.;
  deriv(1,0) = PI*cos(PI*(x[1]+1.)/2.)/2.;
}
*/
/*void Forcing1(TPZVec<REAL> &x, TPZVec<REAL> &disp) {
  //  disp[0] =  -2.;
  disp[0] = 0.;
}
void ExactSolution(TPZVec<REAL> &x, TPZVec<REAL> &u, TPZFMatrix &deriv) {
  u[0] = (x[1]+1.);
  deriv(0,0) = 0.;
  deriv(1,0) = 1.;
}*/

TPZCompMesh *CreateMesh();

int main(){

  int p, h;
  char filename[20];
  char filedx[20];
  cout << "\nOrdem p" << endl;

  cin >> p;

  cout << "\nRefinamento" << endl;

  cin >> h;

  cout << "\nArquivo" << endl;

  cin >> filename;

  cout << "\nArquivoDX" << endl;

  cin >> filedx;

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

  ofstream out(filename), outdx(filedx);

  TPZCompEl::gOrder = 1;
  TPZCompElDisc::gDegree = p;

  gDebug = 0;

  TPZCompMesh *cmesh = CreateMesh();

  TPZGeoMesh *gmesh = cmesh->Reference();


  //  gmesh->Print(out);
  //  cmesh->Print(out);

  TPZAnalysis an(cmesh);
  TPZFStructMatrix full(cmesh);
 
  an.SetStructuralMatrix(full);

  TPZStepSolver step;

  step.SetDirect(ELU);
  an.SetSolver(step);


  /*
  an.Assemble();

  TPZFMatrix mysol(8,1,0.),myres(8,1,0.);
  mysol(1,0) = 1.;
  mysol(3,0) = 0.5;
  mysol(5,0) = 1.;
  mysol(7) = 1.5;
  an.Solver().Matrix()->Multiply(mysol,myres);

  myres.Print("residuo",out);
  */
  an.Run();

  // an.Print( "Nosso primeiro teste",out);

  TPZVec<char *> scalnames(1);
  TPZVec<char *> vecnames(1);
  scalnames[0] = "Solution";
  vecnames[0] = "Derivate";
  an.DefineGraphMesh(2,scalnames,vecnames,filedx);

  an.PostProcess(4);
  an.SetExact(ExactSolution);

  TPZVec<REAL> pos;
  an.PostProcess(pos,out);
  
  delete cmesh;
  delete gmesh;
  return 0;
}

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
  int indices[4][4] = {{0,1,2,3},{0,3,4,5},{0,5,6,7},{0,7,8,1}};
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

    
  for(int i = 0; i < nelem; i++){
    TPZVec<TPZGeoEl *> children, netos, bisnetos, tata1, tata2, tata3;
    cout << "\ngDivide[0] = \n" << gDivide[0];
    if(gDivide[0] == 1)
      {
      elvec[i]->Divide(children);
      cout <<  "\n Primeira divisao \n" ;
      }
    if (gDivide[1] == 1)
      {
      for(int j = 0; j < children.NElements(); j++)
	children[j]->Divide(netos);
      cout <<  "\n Segunda divisao \n" ;
      }
    if(gDivide[2] == 1) 
      {
      for(int k = 0; k < netos.NElements(); k++)
	netos[k]->Divide(bisnetos); 
      cout <<  "\n Terceira divisao \n" ;
      }
    if(gDivide[3] == 1) 
      {
      for(int k = 0; k < bisnetos.NElements(); k++)
	bisnetos[k]->Divide(tata1); 
      cout <<  "\n Quarta divisao \n" ;
      }
    if(gDivide[4] == 1) 
      {
      for(int k = 0; k < tata1.NElements(); k++)
	tata1[k]->Divide(tata2); 
      cout <<  "\n Quinta divisao \n" ;
      }
    if(gDivide[5] == 1) 
      {
      for(int k = 0; k < tata2.NElements(); k++)
	tata2[k]->Divide(tata3); 
      cout <<  "\n Sexta divisao \n" ;
      }

  }  
  
  TPZGeoElBC gbc;

/* Lineare em y:
  TPZGeoElBC gbc1(elvec[0],5,-1,*gmesh); // bottom
  TPZGeoElBC gbc4(elvec[1],6,-2,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-2,*gmesh); // top
  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom
*/

  // bc -1 -> Dirichlet homogeneo
  TPZGeoElBC gbc1(elvec[0],5,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],6,-2,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-2,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-1,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],5,-1,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],6,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],5,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],6,-1,*gmesh); // bottom
  
  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  cmesh->SetDimModel(2);
  
  TPZMaterial *mat;
  mat = new TPZMatPoisson3d(1, 2);
  mat->SetForcingFunction(Forcing1);

  int nstate = 1;
  TPZFMatrix val1(nstate,nstate,0.),val2(nstate,1,0.);
  TPZBndCond *bc[2];
  

  bc[0] = mat->CreateBC(-1,0,val1,val2);
  //  bc[0]->SetForcingFunction(Dirichlet1);

  bc[1] = mat->CreateBC(-2,0,val1,val2);
  //  bc[1]->SetForcingFunction(Dirichlet2);
  
  cmesh->InsertMaterialObject(mat);
  int i;
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);

  TPZGeoElement<TPZShapeCube,TPZGeoCube,TPZRefCube>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeLinear,TPZGeoLinear,TPZRefLinear>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeQuad,TPZGeoQuad,TPZRefQuad>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTriang,TPZGeoTriangle,TPZRefTriangle>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePrism,TPZGeoPrism,TPZRefPrism>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapeTetra,TPZGeoTetrahedra,TPZRefTetrahedra>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  TPZGeoElement<TPZShapePiram,TPZGeoPyramid,TPZRefPyramid>::SetCreateFunction(TPZCompElDisc::CreateDisc);
  //template class TPZGeoElement<TPZShapePoint,TPZGeoPoint,TPZRefPoint>;

  TPZCompElDisc::gInterfaceDimension = 1;
  
  cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}

