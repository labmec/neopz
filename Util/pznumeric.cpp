/** 
 * @file
 * @brief Contains the implementation of the methods to TPZNumeric class.
 */
/***************************************************************************
 pznumeric.cpp  -  description
 -------------------
 begin                : Wed Mar 27 2002
 copyright            : (C) 2002 by Renato Gomes Damas
 email                : rgdamas@fec.unicamp.br
 ***************************************************************************/

/***************************************************************************
 *                                                                         *
 *   This program is free software; you can redistribute it and/or modify  *
 *   it under the terms of the GNU General Public License as published by  *
 *   the Free Software Foundation; either version 2 of the License, or     *
 *   (at your option) any later version.                                   *
 *                                                                         *
 ***************************************************************************/

#include "pznumeric.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>

using namespace std;

REAL TPZNumeric::Norma(const TPZVec<REAL> &vetor){
	const int size = vetor.NElements();
	return sqrt(inner_product(&vetor[0], &vetor[size], &vetor[size], REAL(0.0)));
}
void TPZNumeric::NormalizeVetor3(TPZVec<REAL> &vetor){
	int i;
	REAL norma = TPZNumeric::Norma(vetor);
	for(i=0; i<3; i++) vetor[i] = vetor[i]/norma;
}

TPZNumeric::TPZNumeric(){
}
TPZNumeric::~TPZNumeric(){
}
/** Dada a array[3]armazena sua ordem decrescente, em valor absoluto, em ordem[3]. */
void TPZNumeric::SortArray3(const TPZVec<REAL> &array,int ordem[3]){
    int i;
    REAL vetor[3];
    for(i=0; i<3; i++) vetor[i]=fabs(array[i]);
	
	if (vetor[0] >= vetor[1]) {
		if (vetor[1] >= vetor[2]) {
			ordem[0] = 0;
			ordem[1] = 1;
			ordem[2] = 2;
		}
		else if (vetor[2] >= vetor[0]) {
			ordem[0] = 2;
			ordem[1] = 0;
			ordem[2] = 1;
		}
		else{
			ordem[0] = 0;
			ordem[1] = 2;
			ordem[2] = 1;
		}
	}
	else {
		if (vetor[0] >= vetor[2]) {
			ordem[0] = 1;
			ordem[1] = 0;
			ordem[2] = 2;
		}
		else if (vetor[2] >= vetor[1]) {
			ordem[0] = 2;
			ordem[1] = 1;
			ordem[2] = 0;
		}
		else {
			ordem[0] = 1;
			ordem[1] = 2;
			ordem[2] = 0;
		}
	}
}

void TPZNumeric::SortArray3(TPZVec<REAL> &array){
	int size = array.NElements();
	sort(&array[0], &array[size], greater<REAL>());
}
/** dados dois vetores calcula o produto vetorial. */
void TPZNumeric::ProdVetorial(TPZVec<REAL> &u, TPZVec<REAL> &v, TPZVec<REAL> &result){
	int i;
	REAL aux[3];
	for(i=0; i<3; i++){
		aux[i]= u[(i+1)%3]*v[(i+2)%3] - u[(i+2)%3]*v[(i+1)%3];
	}
	for(i=0; i<3; i++){
		result[i]= aux[i];
	}
}
