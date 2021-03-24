/** 
 * @file
 * @brief Contains the implementation of the methods to TPZNumeric class.
 */

#include "pznumeric.h"
#include <iostream>
#include <cmath>
#include <algorithm>
#include <numeric>
#include <functional>

using namespace std;

template <class Tvar>
Tvar TPZNumeric::Norm(const TPZVec<Tvar> &vetor){
	const int size = vetor.NElements();
    Tvar norma=0;
    for (int i=0; i<size; i++) {
        norma += vetor[i]*vetor[i];
    }
	return sqrt(norma);
}

template <class Tvar>
void TPZNumeric::NormalizeVetor3(TPZVec<Tvar> &vetor){
	int i;
	Tvar norma = TPZNumeric::Norm(vetor);
	for(i=0; i<3; i++) vetor[i] = vetor[i]/norma;
}

template <class Tvar>
void TPZNumeric::NormalizeVetor(TPZVec<Tvar> &vetor){
    int i,n;
    n = vetor.NElements();
    Tvar norma = TPZNumeric::Norm(vetor);
    for(i=0; i< n; i++) vetor[i] = vetor[i]/norma;
}

TPZNumeric::TPZNumeric(){
}
TPZNumeric::~TPZNumeric(){
}

/** Dada a array[3]armazena sua ordem decrescente, em valor absoluto, em ordem[3]. */
template <class Tvar>
void TPZNumeric::SortArray3(const TPZVec<Tvar> &array,int ordem[3]){
    int i;
    Tvar vetor[3];
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

template <class Tvar>
void TPZNumeric::SortArray3(TPZVec<Tvar> &array){
	int size = array.NElements();
	sort(&array[0], &array[size], greater<Tvar>());
}

/** dados dois vetores calcula o produto vetorial. */
template <class Tvar>
void TPZNumeric::ProdVetorial(TPZVec<Tvar> &u, TPZVec<Tvar> &v, TPZVec<Tvar> &result){
	int i;
	Tvar aux[3];
	for(i=0; i<3; i++){
		aux[i]= u[(i+1)%3]*v[(i+2)%3] - u[(i+2)%3]*v[(i+1)%3];
	}
	for(i=0; i<3; i++){
		result[i]= aux[i];
	}
}

template
void TPZNumeric::SortArray3<REAL>(const TPZVec<REAL> &array, int ordem[3]);
template
void TPZNumeric::SortArray3<REAL>(TPZVec<REAL> &array);
template
void TPZNumeric::ProdVetorial<REAL>(TPZVec<REAL> &u, TPZVec<REAL> &v, TPZVec<REAL> &result);
template
void TPZNumeric::NormalizeVetor3<REAL>(TPZVec<REAL> &vetor);
template
void TPZNumeric::NormalizeVetor<REAL>(TPZVec<REAL> &vetor);
template
REAL TPZNumeric::Norm(const TPZVec<REAL> &vetor);


template
void TPZNumeric::SortArray3<Fad<REAL>>(const TPZVec<Fad<REAL>> &array, int ordem[3]);
template
void TPZNumeric::SortArray3<Fad<REAL>>(TPZVec<Fad<REAL>> &array);
template
void TPZNumeric::ProdVetorial<Fad<REAL>>(TPZVec<Fad<REAL>> &u, TPZVec<Fad<REAL>> &v, TPZVec<Fad<REAL>> &result);
template
void TPZNumeric::NormalizeVetor3<Fad<REAL>>(TPZVec<Fad<REAL>> &vetor);
template
void TPZNumeric::NormalizeVetor<Fad<REAL>>(TPZVec<Fad<REAL>> &vetor);
template
Fad<REAL> TPZNumeric::Norm(const TPZVec<Fad<REAL>> &vetor);

