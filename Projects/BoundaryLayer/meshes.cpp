//$Id: meshes.cpp,v 1.5 2006-03-04 15:36:23 tiago Exp $

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

 using namespace pzgeom;
 using namespace pzshape;
 using namespace pzrefine;
 using namespace std;

static REAL epsilon = 1.e-6;;

void OneContinuous(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int continuousindex){

   contset.clear();
   discset.clear();
 
   TPZAdmChunkVector<TPZGeoEl *> &gelvec = gmesh->ElementVec();
   int n = gelvec.NElements();

   contset.insert(gelvec[continuousindex]);
   for(int i = 0; i < n; i++ ){   
      if (i == continuousindex ) continue;
      discset.insert( gelvec[i] );
   } 
 
   TPZManVector<TPZGeoEl*> filhos;
   for(int i = 0; i < h; i++){

      int n = gmesh->NElements();

      for(int j = 0; j < n; j++){
	
	 TPZGeoEl * gel = gmesh->ElementVec()[j];
	 if ( gel->HasSubElement()  ) continue;
	 if ( gel->Dimension() != 2 ) continue;
	 
	 gel->Divide(filhos);

	 if (contset.count(gel)) {
	    int nsons = filhos.NElements();
	    for(int ison = 0; ison < nsons; ison++) 
	       contset.insert(filhos[ison]);
	 }
	
	 if (discset.count(gel)){
	    int nsons = filhos.NElements();
	    for(int ison = 0; ison < nsons; ison++) 
	       discset.insert(filhos[ison]);
	 }
	 
      }
   }


}



void OneDiscontinuous(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int discontinuousindex){

   contset.clear();
   discset.clear();
 
   TPZAdmChunkVector<TPZGeoEl *> &gelvec = gmesh->ElementVec();
   int n = gelvec.NElements();

   discset.insert(gelvec[discontinuousindex]);
   for(int i = 0; i < n; i++ ){   
      if (i == discontinuousindex ) continue;
      contset.insert( gelvec[i] );
   } 
 
   TPZManVector<TPZGeoEl*> filhos;
   for(int i = 0; i < h; i++){

      int n = gmesh->NElements();

      for(int j = 0; j < n; j++){
	
	 TPZGeoEl * gel = gmesh->ElementVec()[j];
	 if ( gel->HasSubElement()  ) continue;
	 if ( gel->Dimension() != 2 ) continue;
	 
	 gel->Divide(filhos);

	 if (contset.count(gel)) {
	    int nsons = filhos.NElements();
	    for(int ison = 0; ison < nsons; ison++) 
	       contset.insert(filhos[ison]);
	 }
	
	 if (discset.count(gel)){
	    int nsons = filhos.NElements();
	    for(int ison = 0; ison < nsons; ison++) 
	       discset.insert(filhos[ison]);
	 }
	 
      }
   }


}


void DiscontinuousOnBoundaryLayer(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h){
   contset.clear();
   discset.clear();

   TPZAdmChunkVector<TPZGeoEl *> &gelvec = gmesh->ElementVec();
   int n = gelvec.NElements();
   TPZManVector<TPZGeoEl*> filhos;
   for(int i = 0; i < h; i++){      
      int n = gmesh->NElements();      
      for(int j = 0; j < n; j++){	
	 TPZGeoEl * gel = gmesh->ElementVec()[j];
	 if ( gel->HasSubElement()  ) continue;
	 if ( gel->Dimension() != 2 ) continue;	 
	 gel->Divide(filhos);
      }
   }

   n = gelvec.NElements();
   for( int i = 0; i < n; i++)
      contset.insert( gmesh->ElementVec()[i] );

   n = gmesh->BCElementVec().NElements();
   TPZGeoEl *bcgel;
   int side;
   TPZStack<TPZGeoElSide> stack;
   for( int i = 0; i < n; i++){
      bcgel = gmesh->BCElementVec()[i].fElement;
      side  = gmesh->BCElementVec()[i].fSide;

      if(bcgel->HasSubElement()) GetSubElements( bcgel, side, stack );
      else{
	 stack.Resize(1);
	 stack[0].SetElement(bcgel);
	 stack[0].SetSide(side);
      }

      int nsubel = stack.NElements();
      for( int is = 0; is < nsubel; is++){
	 TPZGeoEl * disc = stack[is].Element();
	 discset.insert( disc );
	 contset.erase( disc );	 
      }
   }
   

}//end of method

void GetSubElements( TPZGeoEl * gel, int side, TPZStack<TPZGeoElSide> &stack ){

   TPZStack<TPZGeoElSide> subel, subel2;
   gel->GetSubElements2(side, subel);
   int n = subel.NElements();
   for(int i = 0; i < n; i++){
      TPZGeoElSide fine = subel.Pop();
      if ( fine.Element()->HasSubElement() )
      {
	 GetSubElements( fine.Element(), fine.Side(), subel2 );
	 int nsubel2 = subel2.NElements();
	 for(int isubel2 = 0; isubel2 < nsubel2; isubel2++)
	 {
	    stack.Push( subel2.Pop() );
	 }
      }
      else
      {
	 stack.Push(fine);
      }
   }//for
}//end of method



TPZCompMesh * DiscontinuousOnBoundaryLayer(int h){

  REAL co[9][2] = {{0.9,0.9},{0.9,0.},{1.,0.},{1.,0.9},{1.,1.},{0.9,1.},{0.,1.},{0.,0.9},{0.,0.}};
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

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-1,*gmesh); // bottom

//FAZENDO CONTINUO DESCONTINUO

  std::set<TPZGeoEl*> contset, discset;
  TPZVec<TPZGeoEl*> continuous, discontinuous;

  OneContinuous(gmesh, contset, discset, h, 3);//elemento 3 eh continuo - O elemento 3 eh o elemento grande na regiao suave do dominio

  {
     int ncont = contset.size();
     continuous.Resize( ncont );
     continuous.Fill( NULL );
     std::set<TPZGeoEl*>::iterator w, e;
     w = contset.begin();
     e = contset.end();
     for( int i = 0; w != e; w++, i++){
	continuous[i] = *w;
     }
     int ndisc = discset.size();
     discontinuous.Resize( ndisc );
     discontinuous.Fill( NULL );
     w = discset.begin();
     e = discset.end();
     for( int i = 0; w != e; w++, i++){
	discontinuous[i] = *w;
     }
  }

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

  cmesh->AutoBuildContDisc(continuous, discontinuous);
  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}






///////////////////////////////////FUNCTIONS///////////////////////////////////////////

void ForcingFunction(TPZVec<REAL> &pto, TPZVec<REAL> &force){
  REAL x = pto[0];
  REAL y = pto[1];
  
  force.Resize(1);
  REAL num = (1-x)*(1-x)-(1-x)-(1-y)+(1-y)*(1-y);
  REAL den = epsilon*(1-exp(-1/epsilon))*exp((1-x)*(1-y)/epsilon);
  force[0] = 2+(num/den)-x-y;

  //#define LUISFILIPE
#ifdef LUISFILIPE
  REAL antigo = (-(-1.0 + exp(1.0/epsilon)) * epsilon * (-2.0+x+y)
		 +exp( (x+y-x*y)/epsilon ) * (-x + x*x + (-1+y)*y) ) / ( ( -1.0 + exp(1.0/epsilon) ) * epsilon);
//  if ( fabs(force[0] - antigo) > 1.e-10 ){
    cout << "antigo = " << antigo << "  novo = " << force[0] << " x = " << x << " y = " << y << endl;
//  }
#endif

  //*-1 pq o pzpoisson3d implementa -Epsilon Laplac(u) + div(beta*u) = -force  
  force[0] *= -1.0;
}

void ExactSolution(TPZVec<REAL> &pto, TPZVec<REAL> &u, TPZFMatrix &deriv) {

  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  deriv.Resize(2,1);
  
  REAL num = exp(-1./epsilon);
  REAL den = 1. - exp(-1./epsilon);
  u[0] = num/den;
  
  num = -1. * exp( (-1.+x) * (1.-y) / epsilon );
  den = 1. - exp(-1./epsilon);  
  u[0] += num/den;
  
  u[0] += x + y -x*y;
  
  REAL A = epsilon;
  
  num = - exp((-1. + x)*(1. - y)/A);
  den = A*(1. - exp(-1./A));
  
  deriv(0,0) = 1. - num/den;
  deriv(0,0) += - y;
  
  num = exp((-1. + x)*(1. - y)/A) * y;
  den = A * (1. - exp(-1./A));
    
  deriv(0,0) += num/den;
  
  
  deriv(1,0) = 1.;
  num = - exp((-1. + x)*(1. - y)/A);
  den = A * (1. - exp(-1./A));
  deriv(1,0) += num/den;
  deriv(1,0) += - x;
  num = exp((-1 + x)*(1. - y)/A) * x;
  den = A*(1. - exp(-1./A));
  deriv(1,0) += num/den;
  
//  cout << pto << "\t" << u << "\t" << deriv << endl;
  
}

void Dirichlet_X_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u) {
  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  
  REAL num = exp(-1./epsilon);
  REAL den = 1. - exp(-1./epsilon);
  u[0] = num/den;
  
  num = -1. * exp( (-1.+x) * (1.-y) / epsilon );
  den = 1. - exp(-1./epsilon);  
  u[0] += num/den;
  
  u[0] += x + y -x*y;
}

void Dirichlet_Y_IgualA_0(TPZVec<REAL> &pto, TPZVec<REAL> &u) {
  REAL x = pto[0];
  REAL y = pto[1];
  u.Resize(1);
  
  REAL num = exp(-1./epsilon);
  REAL den = 1. - exp(-1./epsilon);
  u[0] = num/den;
  
  num = -1. * exp( (-1.+x) * (1.-y) / epsilon );
  den = 1. - exp(-1./epsilon);  
  u[0] += num/den;
  
  u[0] += x + y -x*y;
}



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

  //Dirichlet nulo
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

  cmesh->AdjustBoundaryElements();
  //  cmesh->CleanUpUnconnectedNodes();
  //  cmesh->ExpandSolution();
  
  return cmesh;
}

#include "TPZInterfaceEl.h"
TPZCompMesh *CreateMeshContDisc(int h) {

  TPZInterfaceElement:: SetCalcStiffContDisc();

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



//FAZENDO CONTINUO DESCONTINUO

  std::set<TPZGeoEl*> contset, discset;
  TPZVec<TPZGeoEl*> continuous, discontinuous;

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-1,*gmesh); // bottom

  OneContinuous(gmesh, contset, discset, h, 3);  
  //OneDiscontinuous(gmesh, contset, discset, h);  

  {
     int ncont = contset.size();
     continuous.Resize( ncont );
     continuous.Fill( NULL );
     std::set<TPZGeoEl*>::iterator w, e;
     w = contset.begin();
     e = contset.end();
     for( int i = 0; w != e; w++, i++){
	continuous[i] = *w;
     }
     int ndisc = discset.size();
     discontinuous.Resize( ndisc );
     discontinuous.Fill( NULL );
     w = discset.begin();
     e = discset.end();
     for( int i = 0; w != e; w++, i++){
	discontinuous[i] = *w;
     }
  }

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

  cmesh->AutoBuildContDisc(continuous, discontinuous);

  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
  
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


/** Malha com refinamento na camada limite 14fev2006 **/
TPZCompMesh * RefinedOnBoundLayer(int h, int ref_uniforme ){

  REAL co[9][2] = {{0.9,0.9},{0.9,0.},{1.,0.},{1.,0.9},{1.,1.},{0.9,1.},{0.,1.},{0.,0.9},{0.,0.}};
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
  
  //criando elementos de contorno da camada limite
  TPZManVector<int,2> nodind(2);
  int index;
  nodind[0] = 2; nodind[1] = 3;
  gmesh->CreateGeoElement(EOned, nodind, -3, index);
  nodind[0] = 3; nodind[1] = 4;
  gmesh->CreateGeoElement(EOned, nodind, -3, index);
  nodind[0] = 4; nodind[1] = 5;
  gmesh->CreateGeoElement(EOned, nodind, -3, index);
  nodind[0] = 5; nodind[1] = 6;
  gmesh->CreateGeoElement(EOned, nodind, -3, index);

  gmesh->BuildConnectivity();
  

  TPZGeoElBC gbc1(elvec[0],4,-1,*gmesh); // bottom
//  TPZGeoElBC gbc2(elvec[0],5,-3,*gmesh); // right
//  TPZGeoElBC gbc3(elvec[1],5,-3,*gmesh); // right
//  TPZGeoElBC gbc4(elvec[1],6,-3,*gmesh); // top
//  TPZGeoElBC gbc5(elvec[2],6,-3,*gmesh); // top
  TPZGeoElBC gbc6(elvec[2],7,-2,*gmesh); // left
  TPZGeoElBC gbc7(elvec[3],7,-2,*gmesh); // left
  TPZGeoElBC gbc8(elvec[3],4,-1,*gmesh); // bottom

//FAZENDO CONTINUO DESCONTINUO

  std::set<TPZGeoEl*> contset, discset;
  TPZVec<TPZGeoEl*> continuous, discontinuous;

  OneContinuousRefinementOnBoundLayer(gmesh, contset, discset, h, ref_uniforme, 3, -3);//elemento 3 eh continuo - O elemento 3 eh o elemento grande na regiao suave do dominio

  {
     int ncont = contset.size();
     continuous.Resize( ncont );
     continuous.Fill( NULL );
     std::set<TPZGeoEl*>::iterator w, e;
     w = contset.begin();
     e = contset.end();
     for( int i = 0; w != e; w++, i++){
	continuous[i] = *w;
     }
     int ndisc = discset.size();
     discontinuous.Resize( ndisc );
     discontinuous.Fill( NULL );
     w = discset.begin();
     e = discset.end();
     for( int i = 0; w != e; w++, i++){
	discontinuous[i] = *w;
     }
  }

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

  cmesh->AutoBuildContDisc(continuous, discontinuous);
//cmesh->AutoBuild();
  cmesh->AdjustBoundaryElements();
  cmesh->CleanUpUnconnectedNodes();
  cmesh->ExpandSolution();
  
  return cmesh;  

}//RefinedOnBoundLayer

/** Malha com refinamento na camada limite 14fev2006 **/
void OneContinuousRefinementOnBoundLayer(TPZGeoMesh *gmesh, std::set<TPZGeoEl*> &contset, std::set<TPZGeoEl*> &discset, int h, int ref_uniforme, int continuousindex, int BLMaterialId){

   contset.clear();
   discset.clear();
 
   TPZAdmChunkVector<TPZGeoEl *> &gelvec = gmesh->ElementVec();
   int n = gelvec.NElements();

   contset.insert(gelvec[continuousindex]);
   for(int i = 0; i < n; i++ ){   
      if (i == continuousindex ) continue;
      discset.insert( gelvec[i] );
   } 
   
  /** Refinamento uniforme */
   TPZManVector<TPZGeoEl*> filhos;
   for(int i = 0; i < ref_uniforme; i++){
     int n = gmesh->NElements();
     for(int j = 0; j < n; j++){
       TPZGeoEl * gel = gmesh->ElementVec()[j];
       if ( gel->HasSubElement()  ) continue;
       if ( gel->Dimension() != 2 ) continue;
       gel->Divide(filhos);
       if (contset.count(gel)) {
         int nsons = filhos.NElements();
         for(int ison = 0; ison < nsons; ison++) 
             contset.insert(filhos[ison]);
         }

       if (discset.count(gel)){
         int nsons = filhos.NElements();
         for(int ison = 0; ison < nsons; ison++) 
             discset.insert(filhos[ison]);
         }       
       }
   }   
 
//   TPZManVector<TPZGeoEl*> filhos;
   for(int i = 0; i < h; i++){

      int n = gmesh->NElements();

      for(int j = 0; j < n; j++){

         TPZGeoEl * gel = gmesh->ElementVec()[j];
         
         if ( gel->HasSubElement()  ) continue;
         if ( gel->Dimension() != 2 ) continue;
         
         int nsides = gel->NSides();
         bool isneighbourtoBL = false;
         for (int is = 0; is < nsides; is++){
           TPZGeoElSide thisside(gel, is);
           TPZGeoElSide neighbour = gel->Neighbour(is);
           if(!neighbour.Exists()) continue;
           while(neighbour != thisside) {
             int matId = neighbour.Element()->MaterialId();
             if (matId == BLMaterialId){
               isneighbourtoBL = true;
               break;
             }
             neighbour = neighbour.Neighbour();
           }//while
           if (isneighbourtoBL) break;
         }//for is         

         if (!isneighbourtoBL) continue;
         
         gel->Divide(filhos);

         if (contset.count(gel)) {
           int nsons = filhos.NElements();
	   for(int ison = 0; ison < nsons; ison++) 
	       contset.insert(filhos[ison]);
	 }
	
	 if (discset.count(gel)){
	    int nsons = filhos.NElements();
	    for(int ison = 0; ison < nsons; ison++) 
	       discset.insert(filhos[ison]);
	 }
	 
      }
   }
}//OneContinuousRefinementOnBoundLayer

