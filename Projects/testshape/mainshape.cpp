#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgengrid.h"

#include "pzl2projection.h"

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"
#include "pzshapecube.h"

#include <iostream>

using namespace pzshape;

void InitialMesh(TPZCompMesh * , int );

int main() {
	ofstream saida("malha.txt");

	// Insere malha geometrica e computacional
	TPZCompMesh *cmesh = 0;
	InitialMesh(cmesh,2);
	cmesh->Print(saida);

	// Funcoes shape
	TPZFMatrix phi(100,1,0.),dphi(1,100,0.);

	TPZVec<int> order(27,3);
	TPZVec<int> nodeids(27,0);
	int i;
	for(i=0; i<27; i++) nodeids[i] = i;

	TPZVec<REAL> point(3,0.2);

	int nc,nshape;

	{
		TPZShapeLinear sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(1,nshape);
		
		sel.Shape (point,nodeids,order,phi,dphi);

		phi.Print("Funcoes de forma uni dimensionais");
		dphi.Print("Derivadas de funcoes de forma uni dimensionais");
	}
	{
		TPZShapeQuad sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(2,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma quadrilaterais");
		dphi.Print("Derivadas de funcoes de forma quadri laterais");
	}
	{
		TPZShapeTriang sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(2,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma triangular");
		dphi.Print("Derivadas de funcoes de forma triangular");
	}
	{
		TPZShapeCube sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(3,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma cubo");
		dphi.Print("Derivadas de funcoes de forma cubo");
	}
	{
		TPZShapePiram sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(3,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma piramide");
		dphi.Print("Derivadas de funcoes de forma piramide");
	}
	{
		TPZShapePrism sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(3,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma prisme");
		dphi.Print("Derivadas de funcoes de forma prisme");
	}
	{
		TPZShapeTetra sel;

		nc = sel.NConnects ();
		nshape = sel.NShapeF (order);

		phi.Resize(nshape,1);
		dphi.Resize(3,nshape);
		
		sel.Shape (point,nodeids, order,phi,dphi);

		phi.Print("Funcoes de forma tetrahedra");
		dphi.Print("Derivadas de funcoes de forma tetrahedra");
	}

	return 0;
}

void InitialMesh(TPZCompMesh *c,int dim) {
	if(dim < 1) {
		cout << "Mesh with bad dimension.";
		return;
	}

	TPZGeoMesh *g = new TPZGeoMesh;
	TPZVec<int> nx(dim);
	TPZVec<REAL> X0(dim), X1(dim);
	for(int i=0;i<dim;i++) {
		nx[0] = 3, nx[1] = 3;
		X0[0] = -1., X0[1] = -1.;
		X1[0] = 1., X1[1] = 1.;
	}
	TPZGenGrid gen(nx,X0,X1,1,0.5);
	gen.Read(*g);
	if(c) delete c;
	c = new TPZCompMesh(g);
	TPZVec<REAL> sol(1,0.);
	TPZAutoPointer<TPZMaterial> material = new TPZL2Projection(1,2,1,sol);
	c->InsertMaterialObject(material);
	TPZCompMesh::SetAllCreateFunctionsDiscontinuous();
	TPZCompEl::SetgOrder(1);
	c->SetDefaultOrder(1);
	c->AutoBuild();
}
