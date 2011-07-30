/** 
 * @file
 * @brief Contains the implementation of the methods to TPZPlane class.
 */
/***************************************************************************
                          pzplane.cpp  -  description
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

#include "pzplane.h"
#include <numeric>
using namespace std;
void MatrixDet(REAL matrix[3][3], REAL &det);
REAL MatrixDet(REAL matrix[3][3]);

TPZPlane::TPZPlane(){
}
TPZPlane::~TPZPlane(){
}
/** Dado três pontos calcula a equação do plano que os contém. */
int TPZPlane::SetPlane(const TPZVec<REAL> &p1, const TPZVec<REAL> &p2,const TPZVec<REAL> &p3){
	REAL matrix[3][3];
	TPZVec<REAL> desloc1(3), desloc2(3);
	TPZVec<REAL> aux(3);
	int i ;
	
	for(i=0; i<3;i++){
		desloc1[i] = p1[i]- p2[i];
		desloc2[i] = p2[i]-p3[i];
	}  
	TPZNumeric::ProdVetorial(desloc1, desloc2, aux);
	REAL norm = inner_product(&aux[0], &aux[3], &aux[0], REAL(0.0));
	if(norm <= 1e-10){
		cerr << "TPZPlane::SetPlane - Erro: Pontos alinhados não é possível determinar um único plano\n";
		return 0;
	}
	else{
		for(i=0; i<3;i++){
			matrix[0][0]=p1[i%3];
			matrix[0][1]=p1[(i+1)%3];
			matrix[1][0]=p2[i%3];
			matrix[1][1]=p2[(i+1)%3];
			matrix[2][0]=p3[i%3];
			matrix[2][1]=p3[(i+1)%3];
			matrix[i][2]=1.0;
			MatrixDet(matrix, aux[i]);
		}
		
		int ordem[3];
		REAL vect[3];
		TPZNumeric::SortArray3(aux,ordem);
		i=ordem[0];
		
		matrix[0][0]=p1[i%3];
		matrix[0][1]=p1[(i+1)%3];
		matrix[1][0]=p2[i%3];
		matrix[1][1]=p2[(i+1)%3];
		matrix[2][0]=p3[i%3];
		matrix[2][1]=p3[(i+1)%3];
		vect[0]=-p1[(i+2)%3];
		vect[1]=-p2[(i+2)%3];
		vect[2]=-p2[(i+2)%3];
		
		REAL matrix2[3][3];
		REAL sol[3];
		
		for(i=0; i<3;i++){
			matrix2[0][i%3]=vect[0];
			matrix2[1][i%3]=vect[1];
			matrix2[2][i%3]=vect[2];
			matrix2[0][(i+1)%3]=matrix[0][(i+1)%3];
			matrix2[1][(i+1)%3]=matrix[1][(i+1)%3];
			matrix2[2][(i+1)%3]=matrix[2][(i+1)%3];
			matrix2[0][(i+2)%3]=matrix[0][(i+2)%3];
			matrix2[1][(i+2)%3]=matrix[1][(i+2)%3];
			matrix2[2][(i+2)%3]=matrix[2][(i+2)%3];
			
			sol[i]=MatrixDet(matrix2)/MatrixDet(matrix);
		}
		fCoef[(ordem[0])%3]=sol[0];
		fCoef[(ordem[0]+1)%3]=sol[1];
		fCoef[(ordem[0]+2)%3]= 1.0;
		fCoef[3]=sol[2];
	}
	/*
	 for (i=0; i<4; i++){
	 cout << fCoef[i] <<"\t";
	 }
	 */
	cout <<endl;
	return 1;	
}
/** Verifica se o ponto[3] pertence ao plano. Se pertencer retorna 1, caso contrário 0.*/
bool TPZPlane::Belongs(const TPZVec<REAL> &ponto){
	REAL aux=0.0;
	int i;
	
	for(i=0;i<3;i++){
		aux = aux + ponto[i]*fCoef[i];
	}
	aux = aux + fCoef[3];
	//cout << "aux= "<< aux <<endl;
	if(fabs(aux)<0.000009){
		//cout << "aux= "<< aux <<endl;
		return true;
	}
	else return false;
}
/** Verifica se o plano coincide com plano formado pelos três pontos passados. Se pertencer retorna 1, caso contrário 0. */
bool TPZPlane::Belongs(const TPZVec<REAL> &ponto1, const TPZVec<REAL> &ponto2, const TPZVec<REAL> &ponto3){
	int aux=0;
	int i;
	aux = aux + Belongs(ponto1);
	aux = aux + Belongs(ponto2);
	aux = aux + Belongs(ponto3);
	TPZVec<REAL> aux1(3), aux2(3);
	if(aux==3){
		/**verificando se os pontos estão alinhados */
		//criando vetores de deslocalmentos entre pontos
		for(i=0; i<3; i++){
			aux1[i]=100.*(ponto1[i]-ponto2[i]);
			aux2[i]=100.*(ponto2[i]-ponto3[i]);
		}
		//verificando se os vetores de deslocamentos tem a mesma direção
		aux=0;    
		TPZNumeric::ProdVetorial(aux1, aux2, aux1);
		for(i=0; i<3; i++){
			if(fabs(aux1[i])<0.0000009) aux++;
		}
		if(aux==3){
			cout << "TPZPlane::Belongs(p1, p2, p3) - Warning: Os 3 pontos de entrada estão alinhados\n";
		}
		return true;
	}
	else return false;
}


/** Calcula o determinante da matriz[3][3]. */
void MatrixDet(REAL matrix[3][3], REAL &det){
	det = MatrixDet(matrix);	
}
/** Calcula o determinante da matriz[3][3]. */
REAL MatrixDet(REAL matrix[3][3]){
	int i;
	REAL aux = 0.0;
	for(i=0; i<3; i++){
		aux = aux +(matrix[i%3][0]*matrix[(i+1)%3][1]*matrix[(i+2)%3][2])-(matrix[(i+2)%3][0]*matrix[(i+1)%3][1]*matrix[i%3][2]);
	}
	return aux;	
}
