//We are use vectors and integration rules routines
#include "pzvec.h"
#include "pzquad.h"

//IO
#include <iostream>
using namespace std;

//Functions declaration
REAL Funcao (TPZVec<REAL> &pt, int i, int j, int k, int p);
REAL Alfa (int i, int j, int k, int p);


int main(){
        int i = 1;
        int j = 1;
        int k = 1;
        int p = 1;
        TPZVec<REAL> point(3,0.);
        REAL weight = 0.;
        REAL integral = 0.;

	//=====1D Rule=====================================
	TPZInt1d ordem1d (p);
        int npoints = ordem1d.NPoints();

        int it;
        for (it=0;it<npoints;it++){
        	ordem1d.Point(it,point,weight);
                integral += weight * Funcao(point,i,j,k,p);
        }
	//=====End of 1D Rule================================

        cout << "Integral = " << integral << endl;
        return 0;
}


//Functions definition
REAL Funcao(TPZVec<REAL> &pt, int i, int j, int k, int p){
	int r,s,t;
	REAL result=0.;
        REAL x = 0.;
        REAL y = 0.;
        REAL z = 0.;

	//To treat 1D and 2D
        int dim = pt.NElements();
        if (dim > 2) z = pt[2];
        else z = 1.;
        if (dim > 1) y = pt[1];
        else y = 1.;
        x = pt[0];

   //Atenção --> A função foi alterada aqui para que 1D e 2D tenham resultados diferentes de zero!
        
	for (r=0;r<=i;r++){
		for(s=0;s<=j;s++){
    	for (t=0;t<=k;t++){
				result +=  Alfa(r,s,t,p) * (pow (x,r) + pow (y,s) + pow (z,t));
   		}
		}
	}
	return result;
};

REAL Alfa(int i, int j, int k, int p){
	int n = i*p*p + j*p + k;
        REAL alfa = sin((REAL)n);
        //cout << "alfa = " << alfa << endl;
	return alfa;
}

