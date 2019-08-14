#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgeoel.h"
#include "pzcompel.h"
#include "pzvec.h"
#include "pzerror.h"
#include <pzmat1dlin.h>
#include <pzmattest3d.h>
#include <pzmattest.h>
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include <pzpoisson3d.h>

//Shape Classes
#include <pzshapepoint.h>
#include <pzshapelinear.h>
#include <pzshapetriang.h>
#include <pzshapequad.h>
#include <pzshapetetra.h>
#include <pzshapepiram.h>
#include <pzshapeprism.h>
#include <pzshapecube.h>

//Geometry Classes
#include <pzgeopoint.h>
#include <TPZGeoLinear.h>
#include <pzgeotriangle.h>
#include <pzgeoquad.h>
#include <pzgeotetrahedra.h>
#include <pzgeopyramid.h>
#include <pzgeoprism.h>
#include <TPZGeoCube.h>

//matrix and analysis classes
#include "pzfmatrix.h"
#include "pzanalysis.h"
#include "pzskylstrmatrix.h"

//solvers
#include "pzsolve.h"
#include "pzstepsolver.h"

//IO
#include <iostream>
using namespace std;
using namespace pzshape;
using namespace pzgeom;
int gDebug = 0;
void ExactSimple3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol);

void CheckRefPattern(string & filename);

int main(){

// 	string filename ("/home/pos/cesar/Piramide_Heman");
// 	CheckRefPattern(filename);
// 	return 0;

  REAL coordinates [4][3]={{0.,0.,0.},{1.,0.,0.},{1.,1.,0.},{0.,1.,0.}};
  TPZVec <REAL> coord (3,0.);
  int i,j,index;
  TPZGeoMesh *gmesh = new TPZGeoMesh();
  for (i=0;i<4;i++){
    for (j=0;j<3;j++){
      coord[j] = coordinates[i][j];

    }
    index = gmesh->NodeVec().AllocateNewElement();
    gmesh->NodeVec()[index] = TPZGeoNode(index,coord,*gmesh);
  }

  TPZVec<int> connect(4,0);
  for (i=0;i<4;i++) connect[i] = i;
  int a = 1;
  TPZGeoEl *gel;
  TPZAutoPointer<TPZRefPattern> unifquad;
  if (a){

    // Crio um elemento sem padrão de refinamento -->> \
    // Está implementado na malha como uniforme (fRefPattern  = 0)
    gel = new TPZGeoElRefPattern <TPZShapeQuad,TPZGeoQuad> (index,connect,1,*gmesh,0/*,unifquad*/);
  }
  else {
    gel = gmesh->CreateGeoElement(EQuadrilateral,connect,1,index);
  }

  //Crio um padrão de refinamento e defino-o como padrão de gel
  unifquad = new TPZRefPattern ("/home/pos/cesar/RefPattern/Quad_Unif.rpt");
  TPZGeoElRefPattern <TPZShapeQuad,TPZGeoQuad> *gelrp = dynamic_cast<TPZGeoElRefPattern<TPZShapeQuad,TPZGeoQuad> *> (gel);
  gelrp->SetRefPattern(unifquad);

  // bc -1 -> Dirichlet at the bottom face of the cube
  TPZGeoElBC gbc1(gel,4,-1,*gmesh);
  //  // bc -2 -> Neumann at the top face of the cube
  TPZGeoElBC gbc2(gel,6,-2,*gmesh);
  gmesh->BuildConnectivity();

  gmesh->Print(cout);
  TPZVec<TPZGeoEl *> subel;
  gel->Divide(subel);

  delete gmesh;

/*  TPZCompMesh *cmesh = new TPZCompMesh(gmesh);
  TPZMatPoisson3d *mat = new TPZMatPoisson3d(1,2);
  TPZVec< REAL > convdir (2,1.);

  mat->SetParameters (1.,0.,convdir);
  TPZFMatrix val1(1,1,0.),val2(1,1,0.);
  TPZBndCond *bc[2];
  bc[0] = mat->CreateBC(-1,0,val1,val2);
  val2(0,0) = 1.;
  bc[1] = mat->CreateBC(-2,1,val1,val2);
  for(i=0; i<2; i++) cmesh->InsertMaterialObject(bc[i]);
  cmesh->InsertMaterialObject(mat);
  cmesh->AutoBuild();

  TPZVec<int> subelindex;
  cmesh->ElementVec()[0]->Divide(0,subelindex);
  ofstream out("malha_quad-rp.out");
  cmesh->Print(out);

  TPZSkylineStructMatrix strskyl(cmesh);
  TPZAnalysis an(cmesh);
  an.SetStructuralMatrix(strskyl);
//  an.SetExact(ExactSimple3D);
  TPZStepSolver *direct = new TPZStepSolver;
  direct->SetDirect(ECholesky);
  an.SetSolver(*direct);

  direct = 0;
  an.Run();
  cmesh->Solution().Print("Vetor Solucao",cout);
  an.Solution().Print("Vetor Solucao",cout);

  delete direct;
  delete cmesh;
  delete unifquad;*/
  return 0;


}

void ExactSimple3D(TPZVec<REAL> &x, TPZVec<REAL> &sol, TPZFMatrix &dsol) {
  sol[0] = x[2];
  dsol(0,0) = 0.;
  dsol(1,0) = 1.;
}


void CheckRefPattern(string & filename){
	TPZRefPattern test (filename);
	test.Mesh()->Print(cout);
}


