/** 
 * @file 
 * @brief Contains the implementation of the methods to TPZLine and TPZFunction classes. 
 */

#include "pzfunction.h"

#include "pznumeric.h"
#include "pzline.h"
#include <iostream>
#include <numeric>
#include <assert.h>

using namespace std;

TPZLine::TPZLine():fPoint(3, 0.0), fDirection(3,1.0){
	fTolerance = 0.0001;
}

TPZLine::~TPZLine(){
}

/** Armazena um ponto da reta e sua direção. */
void TPZLine::SetLine(const TPZVec<REAL> &point, const TPZVec<REAL> &dir){
	assert(point.NElements()==3 && dir.NElements()==3);
	fPoint = point;
	fDirection = dir;
}

/** Verifica se o ponto[3] pertence a reta. */
bool TPZLine::Belongs(const TPZVec<REAL> &point){
	assert(point.NElements()==3);
	int i;
	TPZVec<REAL> vetor(3);
	REAL norma;
	vetor = point;
	//Calculando o vetor deslocamento, do point ao fPoint
	for(i=0; i<3; i++){
		vetor[i]= (fPoint[i] - vetor[i]);
	}
	TPZNumeric::ProdVetorial(vetor,fDirection, vetor);
	norma = inner_product(&vetor[0], &vetor[3], &vetor[0], REAL(0.0));
	if (norma < fTolerance) return true;
	else return false;
}

/** Especifica a tolerância para os cálculos */
void TPZLine::SetTolerance(const REAL &tol){
	fTolerance = tol;
}

/** Fornece a tolerância do cálculo, armazenando em tol. */
REAL TPZLine::GetTolerance(){
	return fTolerance;
}

