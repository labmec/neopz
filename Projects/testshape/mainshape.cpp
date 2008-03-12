#include "pzreal.h"

#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzgengrid.h"
#include "pzintel.h"

#include "pzl2projection.h"

#include "pzshapelinear.h"
#include "pzshapequad.h"
#include "pzshapetriang.h"
#include "pzshapetetra.h"
#include "pzshapepiram.h"
#include "pzshapeprism.h"
#include "pzshapecube.h"

// integration rule
#include <pzquad.h>

#include <iostream>
#include <math.h>

using namespace pzshape;


TPZCompMesh *InitialMesh(int dim,int order,int nsubdiv);
void TestShape(TPZInterpolatedElement *el,int order,ostream &out);
void TestShapeWithPrint(TPZInterpolatedElement *el,int order,ostream &saida);

int main() {
	int i, order;
	int dim = 2, max_order;
	TPZInterpolatedElement *el;
	ofstream saida("malha.txt");

	// Insere malha geometrica e computacional
//	cout << "Ordem de integração = ";
//	cin >> int_order;
	max_order = 32;
	TPZCompMesh *cmesh = InitialMesh(dim,max_order,1);

	// Calcula os valores das funcoes shape no lado e no elemento
	for(i=0;i<cmesh->NElements();i++) {
		el = (TPZInterpolatedElement *)cmesh->ElementVec()[i];
		for(order=0;order<max_order;order++)
			TestShape(el,order,saida);
	}
	saida.close();
	return 0;

	/* Funcoes shape
	TPZFMatrix phi(100,1,0.),dphi(1,100,0.);

	TPZVec<int> order(27,3);
	TPZVec<int> nodeids(27,0);
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
*/
	return 0;
}

TPZCompMesh *InitialMesh(int dim,int order,int nsubdiv) {
	if(dim < 1) {
		cout << "Mesh with bad dimension.";
		return 0;
	}

	TPZGeoMesh *g = new TPZGeoMesh;
	TPZVec<int> nx(dim);
	TPZVec<REAL> X0(dim);
	TPZVec<REAL> X1(dim);
	for(int i=0;i<dim;i++) {
		nx[i] = nsubdiv;
		X0[i] = -3.;
		X1[i] = 3.;
	}
	TPZGenGrid gen(nx,X0,X1,1,0.5);
	gen.Read(*g);
	TPZCompMesh *c = new TPZCompMesh(g);
	TPZVec<REAL> sol(1,0.);
	TPZAutoPointer<TPZMaterial> material = new TPZL2Projection(1,2,1,sol);
	c->InsertMaterialObject(material);
	c->SetDefaultOrder(order);
	c->AutoBuild();
	return c;
}

void TestShapeWithPrint(TPZInterpolatedElement *el,int order,ostream &saida) {
	TPZGeoEl *gel = el->Reference();
	int nsides = gel->NSides();
	int npoints;
	TPZVec<REAL> p(3,0);
	TPZVec<REAL> p_int(3,0);
	TPZFMatrix phi(100,1,0.), dphi(3,100,0.);
	TPZFMatrix phiint(100,1,0.), dphiint(3,100,0.);
	int nsideconnects, nsideshapef;
	int nconnects, nshapef;
	int i, j, k;
	TPZIntPoints *pointIntRule = 0;
	TPZVec<REAL> peso(1000,0);
	TPZTransform transform;

	saida << "Elemento: " << el->Index() << "  Order = " << order << endl;
	nconnects = el->NConnects();
	nshapef = el->NShapeF();
	saida << " " << "connects: " << nconnects << "  n shape: " << nshapef << endl;
	for(i=0;i<nsides;i++) {
		nsideconnects = el->NSideConnects(i);
		nsideshapef = el->NSideShapeF(i);
		saida << "  Side " << i << ": " << "  side connects: " << nsideconnects << "  n side shape: " << nsideshapef << endl;
		transform = gel->SideToSideTransform(i,nsides-1);
		pointIntRule = el->Reference()->CreateSideIntegrationRule(i,order);
		npoints = pointIntRule->NPoints();
		for(j=0;j<npoints;j++) {
			pointIntRule->Point(j,p_int,peso[j]);
			transform.Apply(p_int,p);
			//el->Reference()->X(p,p_int);
			el->SideShapeFunction(i,p_int,phi,dphi);
			el->Shape(p,phiint,dphiint);
			saida << " " << " Point " << j << " (" << p[0] << "," << p[1] << "," << p[2] << ")" << "  peso = " << peso[j] << endl;
			saida << " " << " " << " MATRIX Side Shape: ";
			for(k=0;k<nsideconnects;k++) saida << phi(k,0) << "  ";
			saida << endl;
			saida << " " << " " << " MATRIX Shape:      ";
			for(k=0;k<nconnects;k++) 
				saida << phiint(k,0) << "  ";
			saida << endl;
		}
		saida << endl;
	}
	delete pointIntRule;
}

void TestShape(TPZInterpolatedElement *el,int order,ostream &saida) {
	TPZGeoEl *gel = el->Reference();
	int nsides = gel->NSides();
	int npoints;
	TPZVec<REAL> p(3,0);
	TPZVec<REAL> p_int(3,0);
	TPZFMatrix phi(100,1,0.), dphi(3,100,0.);
	TPZFMatrix phiint(100,1,0.), dphiint(3,100,0.);
	int nsideconnects, nsideshapef;
	int nconnects, nshapef;
	int i, j;
	TPZIntPoints *pointIntRule = 0;
	TPZVec<REAL> peso(1000,0);
	TPZTransform transform;
	bool notequal = false;

	nconnects = el->NConnects();
	nshapef = el->NShapeF();
	TPZVec<int> firstindexinternal(el->NConnects()+1,0);
	for(i=0; i<nsides; i++)
	{
		firstindexinternal[i+1] = firstindexinternal[i]+el->NConnectShapeF(i);
	}
	saida << endl;
	for(i=0;i<nsides;i++) {
		nsideconnects = el->NSideConnects(i);
		nsideshapef = el->NSideShapeF(i);
		TPZVec<int> firstindexside(el->NSideConnects(i)+1,0);
		int isc;
		for(isc=0; isc<nsideconnects; isc++) firstindexside[isc+1] = firstindexside[isc]+el->NConnectShapeF(el->SideConnectIndex(isc,i));
		transform = gel->SideToSideTransform(i,nsides-1);
		pointIntRule = el->Reference()->CreateSideIntegrationRule(i,order);
		npoints = pointIntRule->NPoints();
		for(j=0;j<npoints;j++) {
			pointIntRule->Point(j,p_int,peso[j]);
			transform.Apply(p_int,p);
			el->SideShapeFunction(i,p_int,phi,dphi);
			el->Shape(p,phiint,dphiint);
			// Comparando os valores das funções shape nos correspondentes connects
			int isc;
			for(isc=0; isc<el->NSideConnects(i); isc++)
			{
				int sideconnectindex = el->SideConnectIndex(isc,i);
				for(int il=firstindexside[isc], ilint = firstindexinternal[sideconnectindex];il<firstindexside[isc+1]; il++,ilint++)
				{
					if(fabs(phi(il,0)-phiint(ilint,0)) > 1.e-6)
					{
						notequal = true;
					}
				}
			}
		}
		if(notequal) { 
			saida << "\tValor diferente. Side " << i << endl;
			notequal = false;
		}
	}
	delete pointIntRule;
}
