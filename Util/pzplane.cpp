/** 
 * @file
 * @brief Contains the implementation of the methods to TPZPlane class.
 */

#include "pzplane.h"
#include <numeric>
using namespace std;
void MatrixDet(REAL matrix[3][3], REAL &det);
REAL MatrixDet(REAL matrix[3][3]);

TPZPlane::TPZPlane(){
}
TPZPlane::~TPZPlane(){
}

// Calculate the equation of the plane defined by the given 3 points
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
		cerr << "TPZPlane::SetPlane - Error: Colinear points can't define a single plane\n";
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
	cout <<endl;
	return 1;	
}

// Verify if a point belong to the plane
bool TPZPlane::Belongs(const TPZVec<REAL> &ponto){
	REAL aux=0.0;
	int i;
	
	for(i=0;i<3;i++){
		aux = aux + ponto[i]*fCoef[i];
	}
	aux = aux + fCoef[3];
	if(fabs(aux)<0.000009) {
		return true;
	}
	else return false;
}

// Verify if the plane and the plane formed by the 3 given points are coincident
bool TPZPlane::Belongs(const TPZVec<REAL> &ponto1, const TPZVec<REAL> &ponto2, const TPZVec<REAL> &ponto3){
	int aux=0;
	int i;
	aux = aux + Belongs(ponto1);
	aux = aux + Belongs(ponto2);
	aux = aux + Belongs(ponto3);
	TPZVec<REAL> aux1(3), aux2(3);
    // Verify if the points are colinear
	if(aux==3) {
		for(i=0;i<3;i++) {
			aux1[i]=100.*(ponto1[i]-ponto2[i]);
			aux2[i]=100.*(ponto2[i]-ponto3[i]);
		}
		// Verify if the displacement vectors have the same direction
		aux=0;
		TPZNumeric::ProdVetorial(aux1, aux2, aux1);
		for(i=0;i<3;i++) {
			if(fabs(aux1[i])<0.0000009) aux++;
		}
		if(aux==3) {
			cout << "TPZPlane::Belongs(p1, p2, p3) - Warning: The 3 input points are colinear\n";
		}
		return true;
	}
	else return false;
}

// Calculate the determinant of the matrix
void MatrixDet(REAL matrix[3][3], REAL &det) {
	det = MatrixDet(matrix);	
}

// Calculate the determinant of the matrix
REAL MatrixDet(REAL matrix[3][3]) {
	int i;
	REAL aux = 0.0;
	for(i=0; i<3; i++) {
		aux = aux +(matrix[i%3][0]*matrix[(i+1)%3][1]*matrix[(i+2)%3][2])-(matrix[(i+2)%3][0]*matrix[(i+1)%3][1]*matrix[i%3][2]);
	}
	return aux;	
}
