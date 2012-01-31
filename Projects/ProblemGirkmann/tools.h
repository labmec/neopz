/*
 *
 *Resolver o problema de Girkmann
 *
 *  Created by Agnaldo on 4/29/09.
 */
#ifndef TOOLS_H
#define TOOLS_H

#include "pzvec.h"
#include "pzfmatrix.h"
#include "pzgmesh.h"
#include "pzcmesh.h"
#include "pzcompel.h"

#include "pzgeoquad.h"
#include "TPZGeoLinear.h"
#include "pzgeotriangle.h"
#include "tpzquadrilateral.h"
#include "pzgnode.h"
#include "tpzarc3d.h"
#include "pzgeopoint.h"

#include "tpzcompmeshreferred.h"
#include "TPZRefPattern.h"
#include "tpzgeoelrefpattern.h"
#include "tpzgeoblend.h"

#include "pzelasAXImat.h" 

#include "tpzautopointer.h"
#include "pzbndcond.h"
#include "pzanalysis.h"
#include "TPZParSkylineStructMatrix.h"
#include "pzstepsolver.h"

class tools{
	
public:
	
	tools(REAL Rc, REAL h, REAL alpha, REAL a, REAL b, REAL Rho,REAL young, REAL poisson, bool AnelComPesoLeve);
	
	/**Dados de Entrada do Problema e do materia
	 *Rc (m): Raio da Coroa (Casca)
	 *h (m): Espessura da casca
	 *alpha (Pi rad): Angulo meridional entre a casca e o anel
	 *a (m): Comprimento do anel de enrijecimento
	 *b (m): Largura do anel de enrijecimento
	 *Rho (N/m^3): Peso Especifico do material
	 *young (pa): Modulo de Young
	 *poisson: Coeficinte de Poisson
	 *AnelComPesoLeve: Se igual true => anel tem peso leve. Caso contrario (false) => anel eh casca tem pesos iguais.
	 */
	
	~tools();
	
	/**
	 *Matriz de revolucao
	 *theta (Pi rad): Angulo meridional entre a casca e o anel 
	 */
	TPZFMatrix MatrixR(REAL theta);
	
	/**
	 *Criar malha geometrica
	 *ndiv: Numero de divisões na coroa
	 *ndirectdivL: numero de refinamento direcional a esquerda do elemento da juncao entre a casca e o anel
	 *ndirectdivR: numero de refinamento direcional a direita do elemento da juncao entre a casca e o anel
	 *ndirectdivp: numero de refinamento direcional nos nós da juncao entre a casca e o anel
	 *interface: se interface==true usa-se galerkin descontinuo na interface
	 *RefDirId: refinamento direcional apenas no material de Id = RefDir (RefDir = 1 casca ou 2 anel)
	 Se Refdir = 0 ==> aplicar refinamento direcional nos dois materiais 
	 */
	TPZGeoMesh * MalhaGeoGen(int ndiv, int ndirectdivL,int ndirectdivR,int ndirectdivp, bool interface,  int RefDirId);
	
	/*
	 *metodo para fazer refinamento uniforme
	 */
	void RefinamentoUniforme(TPZGeoMesh &gMesh, int &nh);
	
	/*
	 **AreaBaseAnel (m2): Area da base do anel (cilindro)
	 */
	REAL AreaBaseAnel();
	
	/*
	 **VMaterial (m3): Volume total (casca + anel)
	 *Se o peso do anel não for assumido, entao, VMaterial = Volume da Casca.
	 */
	REAL VMaterial();
	
	/**
	 *Criar malha computacional de material Axis Symmetric
	 *gMesh: malha geometrica
	 *p: ordem dos polinomios
	 *AreaBaseAnel (m2): Area da base do anel (cilindro)
	 */
	TPZCompMesh * MalhaCompGen(TPZGeoMesh * gMesh, int p);
	
	/**
	 *Criar malha computacional com elemnto de interface do material Axis Symmetric
	 *gMesh: malha geometrica
	 *p: ordem dos polinomios
	 *simetric: metodo simetrico (-1.) ou nao simetrico (+1.)
	 *penalidade: termo de penalizacao:  1. (ou outro valor positivo) para metodo penalizado ou 0. metodo nao penalizado
	 */
	TPZCompMesh * MalhaCompMeshWithInterface(TPZGeoMesh * gMesh, int p,REAL simetric, REAL penalidade);
	
	
	/**
     * Resolve o sistema matricial [M]{c}=[b]
	 */
    void SolveSist(TPZAnalysis &an, TPZCompMesh *fCmesh,  int sim);
	
	/**
	 *Pos-processamento, calcular a cortante e o momento 
	 *na junção entre a casca e o anel: M_alpha e Q_alpha
	 do lado do mat1Id (casca) ou mat2Id (anel)
	 *malha:  malha computacional
	 */
	TPZVec<REAL> CalcCortMomento(TPZCompMesh *malha);
	
	/*
	 *Imprime os elementos de interface
	 */
	void PrintInterface(TPZCompMesh  *malha);
	
	/**
	 * Teste os vetores de carga
	 */
	void TesteInterface(TPZCompMesh *cmesh, TPZFMatrix &solution);
	
	/**
	 * Find the indices of the corner connects
	 */
	void CornerConnects(TPZCompMesh *cmesh, std::set<int> &indices, int matid);
		
protected:
	
	/** Atributos: Dados de Entrada do Problema e do materia
	 *fRc (m): Raio da Coroa (Casca)
	 *fh (m): Espessura da casca
	 *falpha (Pi rad): Angulo meridional entre a casca e o anel
	 *fa (m): Comprimento do anel de enrijecimento
	 *fb (m): Largura do anel de enrijecimento
	 *fRho (N/m^3): Peso Especifico do material
	 *fyoung (pa): Modulo de Young
	 *fpoisson: Coeficinte de Poisson
	 *fAnelComPesoLeve: dessidir se o anel tem peso leve ou se o anel eh casca tem pesos iguais 
	 */
	static REAL fRc;
	static REAL fh;
	static REAL falpha;
	static REAL fa;
	static REAL fb;
	static REAL fRho;
	static REAL fyoung;
	static REAL fpoisson;
	static bool fAnelComPesoLeve;
		
};

#endif
