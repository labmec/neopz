/***************************************************************************
                          pzine.cpp  -  description
                             -------------------
    begin                : Tue Apr 16 2002
    copyright            : (C) 2002 by Renato Gomes Damas
    email                : rgdamas@fec.unicamp.br
 ***************************************************************************/
#include "pznumeric.h"
#include "pzline.h"
#include <iostream>
#include <numeric>
#include <assert.h>

TPZLine::TPZLine():fPoint(3, 0.0), fDirection(3,1.0){
  fTolerance = 0.0001;
}
TPZLine::~TPZLine(){
}
/** Armazena um ponto da reta e sua direção. */
void TPZLine::SetLine(const TPZVec<double> &point, const TPZVec<double> &dir){
  assert(point.NElements()==3 && dir.NElements()==3);
  fPoint = point;
  fDirection = dir;
}
/** Verifica se o ponto[3] pertence a reta. */
bool TPZLine::Belongs(const TPZVec<double> &point){
  assert(point.NElements()==3);
  int i;
  TPZVec<double> vetor(3);
  double norma;
  vetor = point;
  //Calculando o vetor deslocamento, do point ao fPoint
  for(i=0; i<3; i++){
    vetor[i]= (fPoint[i] - vetor[i]);
  }
  TPZNumeric::ProdVetorial(vetor,fDirection, vetor);
  norma = inner_product(&vetor[0], &vetor[3], &vetor[0], 0.0);
  if (norma < fTolerance) return true;
  else return false;
}

/** Especifica a tolerância para os cálculos */
void TPZLine::SetTolerance(const double &tol){
	fTolerance = tol;
}
/** Fornece a tolerância do cálculo, armazenando em tol. */
double TPZLine::GetTolerance(){
	return fTolerance;
}
