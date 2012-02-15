
#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzshapequad.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

void ExemploElasticidade(TPZFMatrix &phi, TPZFMatrix &dphi);
void ExemploGenerico2D(TPZFMatrix &phi, TPZFMatrix &dphi);

using namespace pzshape;
using namespace std;

int main() {

	// definir um vetor que contem a ordem de interpolacao para cada lado
	TPZVec<int> order(5,3);

	// definir um vetor com o identificador de cada no de canto
	TPZVec<int> id(4);
	int i;
	for(i=0; i<4; i++) id[i] = i;

	// definir o vetor que contem o ponto onde calcula-se as funcoes de forma
	TPZVec<REAL> pt(2,0.);

	// definir a matriz que ira armazenar as funcoes de forma e derivadas
	TPZFMatrix phi(16,1),dphi(2,16);

	TPZFlopCounter::gCount.Print();
	// calculo das funcoes de forma
	TPZShapeQuad::Shape(pt,id,order,phi,dphi);

	TPZFlopCounter::gCount.Print();
	// impressao das funcoes de forma
	//	phi.Print("Shape function values");
	//	dphi.Print("Shape function derivative values");

	ExemploElasticidade(phi,dphi);

	TPZFlopCounter::gCount.Print();
	ExemploGenerico2D(phi,dphi);

	TPZFlopCounter::gCount.Print();

	return 0;

}

void ExemploElasticidade(TPZFMatrix &phi, TPZFMatrix &dphi) {

	int i;
	// definir um objeto tipo material elastico
	TPZElasticityMaterial matexemplo(1,1.e5,0.3,1.,1.);

	// imprimir o conteudo do material
	//	matexemplo.Print();

	// Preparar os dados para calculo de contribuicao para matriz de rigidez

	// identificar o numero de variaveis de estado
	int nstate = matexemplo.NStateVariables ();

	// identificar o numero de funcoes de forma
	int nshape = phi.Rows();
	
	// um vetor indicando a posicao do ponto
	TPZVec<REAL> x(2,0.);

	// um vetor com o valor da solucao
	TPZVec<REAL> sol(nstate,0.);

	// uma matriz com as derivadas da funcao no ponto
	TPZFMatrix dsol(2,nstate,0.);

	// uma matriz indicando os eixos correspondentes as derivadas da funcao
	TPZFMatrix axes(3,3,0.);
	for(i=0; i<3; i++) axes(i,i)=1.;
	TPZFMatrix jacinv(axes);

	// a matriz de rigidez e a matriz do vetor de carga
	TPZFMatrix ek(nstate*nshape,nstate*nshape,0.),ef(nstate*nshape,1,0.);

        TPZMaterialData data;
        data.x = x;
        data.jacinv = jacinv;
        data.sol[0] = sol;
        data.dsol[0] = dsol;
        data.axes = axes;
        data.phi = phi;
        data.dphix = dphi;

// 	matexemplo.Contribute (x,jacinv,sol,dsol,1.,axes,phi,dphi,ek,ef);
        matexemplo.Contribute (data,1.,ek,ef);

	//	ek.Print("contribuicao para matriz de rigidez");

	//	ef.Print("contribuicao para vetor de carga");

}

void ExemploGenerico2D(TPZFMatrix &phi, TPZFMatrix &dphi) {

	int i;
	// definir um objeto tipo material elastico
	TPZMat2dLin matexemplo(1);

	// inicializar os dados do material
	TPZFMatrix xk(1,1,1),xc(1,2,1.),xf(1,1,1.);

	matexemplo.SetMaterial(xk,xc,xf);

	// a equacao e : u + du/dx + du/dy = 1

	// imprimir o conteudo do material
	//	matexemplo.Print();

	// Preparar os dados para calculo de contribuicao para matriz de rigidez

	// identificar o numero de variaveis de estado
	int nstate = matexemplo.NStateVariables ();

	// identificar o numero de funcoes de forma
	int nshape = phi.Rows();
	
	// um vetor indicando a posicao do ponto
	TPZVec<REAL> x(2,0.);

	// um vetor com o valor da solucao
	TPZVec<REAL> sol(nstate,0.);

	// uma matriz com as derivadas da funcao no ponto
	TPZFMatrix dsol(2,nstate,0.);

	// uma matriz indicando os eixos correspondentes as derivadas da funcao
	TPZFMatrix axes(3,3,0.);
	for(i=0; i<3; i++) axes(i,i)=1.;
	TPZFMatrix jacinv(axes);

	// a matriz de rigidez e a matriz do vetor de carga
	TPZFMatrix ek(nstate*nshape,nstate*nshape,0.),ef(nstate*nshape,1,0.);

        TPZMaterialData data;
        data.x = x;
        data.jacinv = jacinv;
        data.sol[0] = sol;
        data.dsol[0] = dsol;
        data.axes = axes;
        data.phi = phi;
        data.dphix = dphi;

// 	matexemplo.Contribute (x,jacinv,sol,dsol,1.,axes,phi,dphi,ek,ef);
        matexemplo.Contribute (data,1.,ek,ef);

	//	ek.Print("contribuicao para matriz de rigidez");

	//	ef.Print("contribuicao para vetor de carga");

}
