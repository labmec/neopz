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
void TestShapeIsLinear(TPZInterpolatedElement *el,int order,ostream &saida);

int main() {
	int i, order;
	int dim = 2, max_order;
	TPZInterpolatedElement *el;
	ofstream saida("malha.txt");

	// Insere malha geometrica e computacional
//	cout << "Ordem de integração = ";
//	cin >> int_order;
	max_order = 15;
	TPZCompMesh *cmesh = InitialMesh(dim,max_order,1);

	// Calcula os valores das funcoes shape no lado e no elemento
	for(i=0;i<cmesh->NElements();i++) {
		el = (TPZInterpolatedElement *)cmesh->ElementVec()[i];
		for(order=0;order<max_order;order++)
			TestShapeIsLinear(el,order,saida);
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
	int side, con, point, shape, conside;
	TPZIntPoints *pointIntRule = 0;
	TPZVec<REAL> peso(1000,0);
	TPZTransform transform;

	nconnects = el->NConnects();
	nshapef = el->NShapeF();
	TPZVec<int> firstindexinternal(nconnects+1,0);
	for(con=0;con<nconnects;con++) 
		firstindexinternal[con+1] = firstindexinternal[con]+el->NConnectShapeF(con);
	// Imprime dados do elemento, numero de connects e funçoes shape dependendo da ordem escolhida
	saida << "Elemento: " << el->Index() << "  Order = " << order << endl;
	saida << " " << "connects: " << nconnects << "  n shape: " << nshapef << endl;
	for(side=0;side<nsides;side++) {
		TPZVec<int> firstindexside(nsideconnects+1,0);
		nsideconnects = el->NSideConnects(side);
		nsideshapef = el->NSideShapeF(side);
		for(conside=0;conside<nsideconnects;conside++) {
			firstindexside[conside+1] = firstindexside[conside] + el->NConnectShapeF(el->SideConnectIndex(conside,side));
		}
		saida << "  Side " << side << ": " << "  side connects: " << nsideconnects << "  n side shape: " << nsideshapef << endl;
		transform = gel->SideToSideTransform(side,(nsides-1));
		pointIntRule = el->Reference()->CreateSideIntegrationRule(side,order);
		npoints = pointIntRule->NPoints();
		for(point=0;point<npoints;point++) {
			pointIntRule->Point(point,p_int,peso[point]);
			transform.Apply(p_int,p);
			el->SideShapeFunction(side,p_int,phi,dphi);
			el->Shape(p,phiint,dphiint);
			saida << " " << " Point " << point << " (" << p[0] << "," << p[1] << "," << p[2] << ")" << "  peso = " << peso[point] << endl;
			saida << " " << " " << " MATRIX Side Shape: ";
			for(conside=0;conside<nsideconnects;conside++) {
				nshapef = el->NConnectShapeF(el->SideConnectIndex(conside,side));
				for(shape=0;shape<nshapef;shape++) 
					saida << phi((firstindexside[conside]+shape),0) << "  ";
			}
			saida << endl << " " << " " << " MATRIX Shape:      ";
			for(conside=0;conside<nsideconnects;conside++) {
				nshapef = el->NConnectShapeF(el->SideConnectIndex(conside,side));
				for(shape=0;shape<nshapef;shape++) 
					saida << phiint((firstindexinternal[el->SideConnectIndex(conside,side)]+shape),0) << "  ";
			}
			saida << endl;
			// Repassando e verificando se as diferencas sao nao nulas
			for(conside=0;conside<nsideconnects;conside++) {
				nshapef = el->NConnectShapeF(el->SideConnectIndex(conside,side));
				for(shape=0;shape<nshapef;shape++) {
					if(!IsZero(phi((firstindexside[conside]+shape),0)-phiint((firstindexinternal[el->SideConnectIndex(conside,side)]+shape),0)))
						saida << "  Ops -> " << side << " " << conside << " " << shape << " side,connectside,shapefconnectside." << endl;
				}
			}
		}
		saida << endl;
	}
	delete pointIntRule;
}

void TestShapeIsLinear(TPZInterpolatedElement *el,int order,ostream &saida) {
	TPZGeoEl *gel = el->Reference();
	int nsides = gel->NSides();
	int npoints;
	TPZVec<REAL> p1(3,0.);
	TPZVec<REAL> p2(3,0.);
	TPZVec<REAL> p3(3,0.);
	TPZVec<REAL> p_int(3,0.);
	TPZFMatrix phi1(100,1,0.), dphi1(3,100,0.);
	TPZFMatrix phi2(100,1,0.), dphi2(3,100,0.);
	TPZFMatrix phi3(100,1,0.), dphi3(3,100,0.);
	int nsideconnects, nsideshapef;
	int nconnects, nshapef;
	int side, con, point, shape, conside;
	TPZIntPoints *pointIntRule = 0;
	TPZVec<REAL> peso(1000,0.);
	TPZTransform transform;

	nconnects = el->NConnects();
	nshapef = el->NShapeF();
	TPZVec<int> firstindexinternal(nconnects+1,0);
	for(con=0;con<nconnects;con++) 
		firstindexinternal[con+1] = firstindexinternal[con]+el->NConnectShapeF(con);
	// Imprime dados do elemento, numero de connects e funçoes shape dependendo da ordem escolhida
	saida << "Elemento: " << el->Index() << "  Order = " << order << endl;
	saida << " " << "connects: " << nconnects << "  n shape: " << nshapef << endl;
	for(side=el->NCornerConnects();side<nsides-1;side++) {
		TPZVec<REAL> n(3,0.);
		int i;
		REAL fator = -1.;
		if(!(side%2)) { n[1] = 1.;
			for(i=0;i<side/2;i++)
				n[1] *= fator;
		}
		else { n[0] = 1.;
			for(i=0;i<(side/2 - 1);i++)
				n[1] *= fator;
		}
		TPZVec<int> firstindexside(nsideconnects+1,0);
		nsideconnects = el->NSideConnects(side);
		nsideshapef = el->NSideShapeF(side);
//		for(conside=0;conside<nsideconnects;conside++) {
//			firstindexside[conside+1] = firstindexside[conside] + el->NConnectShapeF(el->SideConnectIndex(conside,side));
//		}
		saida << "  Side " << side << ": " << "  side connects: " << nsideconnects << "  n side shape: " << nsideshapef << endl;
		transform = gel->SideToSideTransform(side,(nsides-1));
		pointIntRule = el->Reference()->CreateSideIntegrationRule(side,order);
		npoints = pointIntRule->NPoints();
		for(point=0;point<npoints;point++) {
			pointIntRule->Point(point,p_int,peso[point]);
			transform.Apply(p_int,p1);
			for(i=0;i<p2.NElements();i++) p2[i] = p1[i]+n[i];
			for(i=0;i<p3.NElements();i++) p3[i] = p2[i]+n[i];
			el->Shape(p1,phi1,dphi1);
			el->Shape(p2,phi2,dphi2);
			el->Shape(p3,phi3,dphi3);
			saida << " " << " Point 1 " << point << " (" << p1[0] << "," << p1[1] << "," << p1[2] << ")" << "  peso = " << peso[point] << endl;
			saida << " " << " Point 2 " << point << " (" << p2[0] << "," << p2[1] << "," << p2[2] << ")" << endl;
			saida << " " << " Point 3 " << point << " (" << p3[0] << "," << p3[1] << "," << p3[2] << ")" << endl;
//			saida << " " << " " << " MATRIX Side Shape: ";
//			for(conside=0;conside<nsideconnects;conside++) {
//				nshapef = el->NConnectShapeF(el->SideConnectIndex(conside,side));
//				for(shape=0;shape<nshapef;shape++) 
//					saida << phi1((firstindexside[conside]+shape),0) << "  ";
//			}
//			saida << endl << " " << " " << " MATRIX Shape:      ";
//			for(conside=0;conside<nsideconnects;conside++) {
//				nshapef = el->NConnectShapeF(el->SideConnectIndex(conside,side));
//				for(shape=0;shape<nshapef;shape++) 
//					saida << phi2((firstindexinternal[el->SideConnectIndex(conside,side)]+shape),0) << "  ";
//			}
			saida << endl;
		}
		saida << endl;
	}
	delete pointIntRule;
}

// Philippe 
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
		for(isc=0; isc<nsideconnects; isc++)
			firstindexside[isc+1] = firstindexside[isc]+el->NConnectShapeF(el->SideConnectIndex(isc,i));
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
			for(isc=0; isc<nsideconnects; isc++)
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
