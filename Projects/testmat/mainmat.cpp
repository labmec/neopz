
#include "pzelasmat.h"
#include "pzmat2dlin.h"
#include "pzshapequad.h"
#include "pzfmatrix.h"
#include "pzvec.h"
#include "pzreal.h"

void ExemploElasticidade(TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);
void ExemploGenerico2D(TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi);

using namespace pzshape;
using namespace std;

int main() {
    
    // Define a vector with polynomial order associated with each side
    TPZVec<int> order(5,3);
    
    // define an id value for each vertex
    TPZVec<int64_t> id(4);
    int i;
    for(i=0; i<4; i++) id[i] = i;
    
    // define the point in parameter space
    TPZVec<REAL> pt(2,0.);
    
    // matrices which will contain the shape functions and their derivatives
    TPZFMatrix<REAL> phi(16,1),dphi(2,16);
    
    // Computing the shape functions
    TPZShapeQuad::Shape(pt,id,order,phi,dphi);
    
    // impressao das funcoes de forma
    phi.Print("Shape function values");
    dphi.Print("Shape function derivative values");
    
    // using the shape functions for an elasticity problem
    ExemploElasticidade(phi,dphi);
    
    // using the shape functions for a material which represents a generic 2d linear system of pdes
    ExemploGenerico2D(phi,dphi);
    
    return 0;
    
}

void ExemploElasticidade(TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    
    int i;
    // Create the elasticity object
    TPZElasticityMaterial matexemplo(1,1.e5,0.3,1.,1.);
    
    // Print the material object
    matexemplo.Print();
    
    // Preparar os dados para calculo de contribuicao para matriz de rigidez
    
    // identify the number of state variables of the material (2)
    int nstate = matexemplo.NStateVariables ();
    
    // identify the number of shape functions
    int nshape = phi.Rows();
    
    // position of the integration point
    TPZVec<REAL> x(2,0.);
    
    // um vetor com o valor da solucao
    TPZManVector<STATE> sol(nstate,(STATE)(0.));
    // uma matriz com as derivadas da funcao no ponto
    TPZFNMatrix<15,STATE> dsolu(2,nstate);
    dsolu.Zero();
    
    // uma matriz indicando os eixos correspondentes as derivadas da funcao
    TPZFMatrix<REAL> axes(3,3,0.);
    for(i=0; i<3; i++) axes(i,i)=1.;
    TPZFMatrix<REAL> jacinv(axes);
    
    // a matriz de rigidez e a matriz do vetor de carga
    TPZFMatrix<STATE> ek(nstate*nshape,nstate*nshape,0.), ef(nstate*nshape,1,0.);
    
    TPZMaterialData data;
    data.x = x;
    data.jacinv = jacinv;
    data.sol[0] = sol;
    data.dsol[0] = dsolu;
    data.axes = axes;
    data.phi = phi;
    data.dphix = dphi;
    
    // 	matexemplo.Contribute (x,jacinv,sol,dsol,1.,axes,phi,dphi,ek,ef);
    matexemplo.Contribute (data,1.,ek,ef);
    
    ek.Print("Contribution to the stiffness matrix");
    
    ef.Print("Contribution to the right hand side");
    
}

void ExemploGenerico2D(TPZFMatrix<REAL> &phi, TPZFMatrix<REAL> &dphi) {
    
    int i;
    // Object of TPZMat2dLin
    TPZMat2dLin matexemplo(1);
    
    // Initialize the data structure
    TPZFMatrix<STATE> xkk(1,1,1),xc(1,2,1.),xf(1,1,1.);
    
    matexemplo.SetMaterial(xkk,xc,xf);
    
    // the equation is : u + du/dx + du/dy = 1
    
    // imprimir o conteudo do material
    matexemplo.Print();
    
    // Prepare the data structure to call Contribute
    
    // identificar o numero de variaveis de estado
    int nstate = matexemplo.NStateVariables ();
    
    // identificar o numero de funcoes de forma
    int nshape = phi.Rows();
    
    // um vetor indicando a posicao do ponto
    TPZVec<REAL> x(2,0.);
    
    // um vetor com o valor da solucao
    TPZVec<STATE> sol(nstate,0.);
    
    // uma matriz com as derivadas da funcao no ponto
    TPZFMatrix<STATE> dsol(2,nstate,0.);
    
    // uma matriz indicando os eixos correspondentes as derivadas da funcao
    TPZFMatrix<REAL> axes(3,3,0.);
    for(i=0; i<3; i++) axes(i,i)=1.;
    TPZFMatrix<REAL> jacinv(axes);
    
    // a matriz de rigidez e a matriz do vetor de carga
    TPZFMatrix<STATE> ek(nstate*nshape,nstate*nshape,0.),ef(nstate*nshape,1,0.);
    
    TPZMaterialData data;
    data.x = x;
    data.jacinv = jacinv;
    data.sol[0] = sol;
    data.dsol[0] = dsol;
    data.axes = axes;
    data.phi = phi;
    data.dphix = dphi;
    
    matexemplo.Contribute (data,1.,ek,ef);
    
    ek.Print("contribuicao para matriz de rigidez");
    
    ef.Print("contribuicao para vetor de carga");
    
}
