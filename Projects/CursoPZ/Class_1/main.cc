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

	int i = 2;
  int j = 2;
  int k = 2;
  int p = 6;
  TPZVec<REAL> point(1,0.);
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
  cout << "1D - Integral = " << integral << endl;
	
	//=====2D Triangle Rule==============================
	point.Resize(2);
	TPZIntTriang ordem2dt (p);
  npoints = ordem2dt.NPoints();
  for (it=0;it<npoints;it++){
   	ordem2dt.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 2D Triangle Rule=======================
  cout << "2D - Triangle Integral = " << integral << endl;

	//=====2D Quad Rule==================================
	TPZIntQuad ordem2dq (p,p);
  npoints = ordem2dq.NPoints();
  for (it=0;it<npoints;it++){
   	ordem2dq.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 2D Quad Rule===========================
  cout << "2D - Quad Integral = " << integral << endl;

	//=====3D Tetra Rule==================================
	point.Resize(3);
	TPZIntTetra3D ordem3dt (p);
  npoints = ordem3dt.NPoints();
  for (it=0;it<npoints;it++){
   	ordem3dt.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 3D Tetra Rule==========================
  cout << "3D - Tetra Integral = " << integral << endl;

	//=====3D Pyramid Rule==================================
	TPZIntPyram3D ordem3dpy (p);
  npoints = ordem3dpy.NPoints();
  for (it=0;it<npoints;it++){
   	ordem3dpy.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 3D Pyramid Rule==========================
  cout << "3D - Pyramid Integral = " << integral << endl;
  
	//=====3D Prism Rule===================================
	TPZIntPrism3D ordem3dpr (p);
  npoints = ordem3dpr.NPoints();
  for (it=0;it<npoints;it++){
   	ordem3dpr.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 3D Prism Rule============================
  cout << "3D - Prism Integral = " << integral << endl;

	//=====3D Hexa Rule====================================
	TPZIntCube3D ordem3dc (p);
  npoints = ordem3dc.NPoints();
  for (it=0;it<npoints;it++){
   	ordem3dc.Point(it,point,weight);
    integral += weight * Funcao(point,i,j,k,p);
  }
	//=====End of 3D Hexa Rule=============================
  cout << "3D - Hexa Integral = " << integral << endl;
  
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

  //cout << "Point coordinate:\n";
	//cout << "x = " << x << "  y = " << y << "  z = " << z << endl;
  
	for (r=0;r<=i;r++){
		for(s=0;s<=j;s++){
    	for (t=0;t<=k;t++){
				result +=  Alfa(r,s,t,p) * (pow (x,r) * pow (y,s) * pow (z,t));
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

