/**
 * \file
 * @brief Contains implementations of the TPZMatHyperElastic methods.
 */
#include "pzmathyperelastic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
//#include "pztempmat.h"
#include "pzerror.h"
#include <math.h>

#include <cmath>
using namespace std;

TPZMatHyperElastic::TPZMatHyperElastic(int nummat,REAL e,REAL mu,REAL nu,
									   REAL lambda,REAL coef1,REAL coef2,REAL coef3) : TPZMaterial(nummat)
{
	
	fE = e;
	fMu = mu;
	fNu = .5*fE/(1.+fMu);
	fLambda = 2.*fNu*fMu/(1.-2.*fMu);
	if(coef1 != -1.0) fCoef1 = coef1;
	else fCoef1 = .25*fLambda;//=a
	if(coef2 != -1.0) fCoef2 = coef2;
	else fCoef2 = -(.5*fLambda+fNu);//=b
	if(coef3 != -1.0) fCoef3 = coef3;
	else fCoef3 = .5*fNu;//=c
	fE1[0] = 2.;
	fE1[1] = 2.;
	fE1[2] = 2.;
	fE5[0] = 2.;
	fE5[1] = 2.;
	fE5[2] = 2.;
	fE9[0] = 2.;
	fE9[1] = 2.;
	fE9[2] = 2.;
}

TPZMatHyperElastic::~TPZMatHyperElastic() {
}

int TPZMatHyperElastic::NStateVariables() {
	return 3;//3 deslocamentos
}

void TPZMatHyperElastic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	TPZMaterial::Print(out);
}

void TPZMatHyperElastic::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<REAL> &ek,TPZFMatrix<REAL> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	// TPZFMatrix<REAL> &dphiL = data.dphixl;
	// TPZFMatrix<REAL> &dphiR = data.dphixr;
	TPZFMatrix<REAL> &phi = data.phi;
	// TPZFMatrix<REAL> &phiL = data.phil;
	// TPZFMatrix<REAL> &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZFMatrix<REAL> &dsol=data.dsol[0];
	// TPZFMatrix<REAL> &dsolL=data.dsoll;
	// TPZFMatrix<REAL> &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix<REAL> &daxesdksi=data.daxesdksi;
	// TPZFMatrix<REAL> &axes=data.axes;
	
	if(fForcingFunction) {
		TPZManVector<REAL> res(3);
		fForcingFunction->Execute(x,res);
		fXf[0] = res[0];
		fXf[1] = res[1];
		fXf[2] = res[2];
	}
#ifdef _AUTODIFF
	TFad<9, TFad<9,REAL> > U;
	ComputeEnergy(fLambda,fNu,dsol,U);
	//   U2 = U*TFad<9,REAL>(weight);
	//   cout <<  "Before multiplication\n" << U << endl;
	//   cout <<  "weight  " << weight <<  "\nAfter multiplication\n" << U2 << endl;
	int nshape = phi.Rows();
	int ish,jsh, i,j;
	for(ish=0; ish<nshape; ish++) {
		for(i=0; i<3; i++) {
			ef(ish*3+i) -= (U.val().d(i)*dphi(0,ish)+U.val().d(3+i)*dphi(1,ish)+U.val().d(6+i)*dphi(2,ish))*weight;
			for(jsh=0; jsh<nshape; jsh++) {
				for(j=0; j<3; j++) {
					ek(ish*3+i,jsh*3+j) += (
											U.d(i).d(j  )*dphi(0,ish)*dphi(0,jsh)+U.d(i+3).d(j  )*dphi(1,ish)*dphi(0,jsh)+U.d(i+6).d(j  )*dphi(2,ish)*dphi(0,jsh)+
											U.d(i).d(j+3)*dphi(0,ish)*dphi(1,jsh)+U.d(i+3).d(j+3)*dphi(1,ish)*dphi(1,jsh)+U.d(i+6).d(j+3)*dphi(2,ish)*dphi(1,jsh)+
											U.d(i).d(j+6)*dphi(0,ish)*dphi(2,jsh)+U.d(i+3).d(j+6)*dphi(1,ish)*dphi(2,jsh)+U.d(i+6).d(j+6)*dphi(2,ish)*dphi(2,jsh)
											)*weight;
				}
			}
		}
	}
	
#else
	int i;
	REAL global[3][3][9];
	REAL ux,uy,uz,vx,vy,vz,wx,wy,wz;
	ux = dsol(0,0);
	uy = dsol(1,0);
	uz = dsol(2,0);
	
	vx = dsol(0,1);
	vy = dsol(1,1);
	vz = dsol(2,1);
	
	wx = dsol(0,2);
	wy = dsol(1,2);
	wz = dsol(2,2);
	
	for(i=0; i<3; i++) fK2[i][i] = 0.;
	fK2[0][0] = 0.;
	fK2[0][1] =  wz+1.;
	fK2[0][2] = -vz;
	fK2[1][0] = -wz-1.;
	fK2[1][1] = 0.;
	fK2[1][2] =  uz;
	fK2[2][0] =  vz;
	fK2[2][1] = -uz;
	fK2[2][2] = 0.;
	
	for(i=0; i<3; i++) fK3[i][i] = 0.;
	fK3[0][1] = -wy;
	fK3[0][2] =  vy+1.;
	fK3[1][0] =  wy;
	fK3[1][2] = -uy;
	fK3[2][0] = -vy-1.;
	fK3[2][1] =  uy;
	
	for(i=0; i<3; i++) fK4[i][i] = 0.;
	fK4[0][1] = -wz-1.;
	fK4[0][2] =  vz;
	fK4[1][0] =  wz+1.;
	fK4[1][2] = -uz;
	fK4[2][0] = -vz;
	fK4[2][1] =  uz;
	
	for(i=0; i<3; i++) fK6[i][i] = 0.;
	fK6[0][1] =  wx;
	fK6[0][2] = -vx;
	fK6[1][0] = -wx;
	fK6[1][2] =  ux+1.;
	fK6[2][0] =  vx;
	fK6[2][1] = -ux-1.;
	
	for(i=0; i<3; i++) fK7[i][i] = 0.;
	fK7[0][1] =  wy;
	fK7[0][2] = -vy-1.;
	fK7[1][0] = -wy;
	fK7[1][2] =  uy;
	fK7[2][0] =  vy+1.;
	fK7[2][1] = -uy;
	
	for(i=0; i<3; i++) fK8[i][i] = 0.;
	fK8[0][1] = -wx;
	fK8[0][2] =  vx;
	fK8[1][0] =  wx;
	fK8[1][2] = -ux-1.;
	fK8[2][0] = -vx;
	fK8[2][1] =  ux+1.;
	
	fGradDetF[0][0] = (vy+1.)*(wz+1.) - wy*vz;//dux
	fGradDetF[1][0] = wx*vz - vx*(wz+1.);     //duy
	fGradDetF[2][0] = vx*wy - wx*(vy+1.);     //duz
	fGradDetF[0][1] = wy*uz - uy*(wz+1.);     //dvx
	fGradDetF[1][1] = (ux+1.)*(wz+1.) - wx*uz;//dvy
	fGradDetF[2][1] = wx*uy - (ux+1.)*wy;     //dvz
	fGradDetF[0][2] = uy*vz - (vy+1.)*uz;     //dwx
	fGradDetF[1][2] = vx*uz - (ux+1.)*vz;     //dwy
	fGradDetF[2][2] = (ux+1.)*(vy+1.) - vx*uy;//dwz
	
	fGradtrC[0][0] = 2.*(ux+1.);//dux
	fGradtrC[1][0] = 2.*uy;     //duy     | dux dvx dwx |
	fGradtrC[2][0] = 2.*uz;     //duz   : | duy dvy dwy |
	fGradtrC[0][1] = 2.*vx;     //dvx     | duz dvz dwz |
	fGradtrC[1][1] = 2.*(vy+1.);//dvy
	fGradtrC[2][1] = 2.*vz;     //dvz
	fGradtrC[0][2] = 2.*wx;     //dwx
	fGradtrC[1][2] = 2.*wy;     //dwy
	fGradtrC[2][2] = 2.*(wz+1.);//dwz
	
	//foi transposta da original
	fL1[0][0] = fGradDetF[0][0]*fGradDetF[0][0];//dux*dux
	fL1[0][1] = fGradDetF[0][0]*fGradDetF[0][1];//dux*dvx
	fL1[0][2] = fGradDetF[0][0]*fGradDetF[0][2];//dux*dwx
	fL1[1][0] = fGradDetF[0][1]*fGradDetF[0][0];//dvx*dux
	fL1[1][1] = fGradDetF[0][1]*fGradDetF[0][1];//dvx*dvx
	fL1[1][2] = fGradDetF[0][1]*fGradDetF[0][2];//dvx*dwx
	fL1[2][0] = fGradDetF[0][2]*fGradDetF[0][0];//dwx*dux
	fL1[2][1] = fGradDetF[0][2]*fGradDetF[0][1];//dwx*dvx
	fL1[2][2] = fGradDetF[0][2]*fGradDetF[0][2];//dwx*dwx
	
	fL2[0][0] = fGradDetF[0][0]*fGradDetF[1][0];//duy*dux
	fL2[0][1] = fGradDetF[0][0]*fGradDetF[1][1];//dvy*dux
	fL2[0][2] = fGradDetF[0][0]*fGradDetF[1][2];//dwy*dux
	fL2[1][0] = fGradDetF[0][1]*fGradDetF[1][0];//duy*dvx
	fL2[1][1] = fGradDetF[0][1]*fGradDetF[1][1];//dvy*dvx
	fL2[1][2] = fGradDetF[0][1]*fGradDetF[1][2];//dwy*dvx
	fL2[2][0] = fGradDetF[0][2]*fGradDetF[1][0];//duy*dwx
	fL2[2][1] = fGradDetF[0][2]*fGradDetF[1][1];//dvy*dwx
	fL2[2][2] = fGradDetF[0][2]*fGradDetF[1][2];//dwy*dwx
	
	fL3[0][0] = fGradDetF[0][0]*fGradDetF[2][0];//duz*dux
	fL3[0][1] = fGradDetF[0][0]*fGradDetF[2][1];//dvz*dux
	fL3[0][2] = fGradDetF[0][0]*fGradDetF[2][2];//dwz*dux
	fL3[1][0] = fGradDetF[0][1]*fGradDetF[2][0];//duz*dvx
	fL3[1][1] = fGradDetF[0][1]*fGradDetF[2][1];//dvz*dvx
	fL3[1][2] = fGradDetF[0][1]*fGradDetF[2][2];//duz*dwx
	fL3[2][0] = fGradDetF[0][2]*fGradDetF[2][0];//dvz*dwx
	fL3[2][1] = fGradDetF[0][2]*fGradDetF[2][1];//dwz*dvx
	fL3[2][2] = fGradDetF[0][2]*fGradDetF[2][2];//dwz*dwx
	
	fL4[0][0] = fGradDetF[1][0]*fGradDetF[0][0];//dux*duy
	fL4[0][1] = fGradDetF[1][0]*fGradDetF[0][1];//dvx*duy
	fL4[0][2] = fGradDetF[1][0]*fGradDetF[0][2];//dwx*duy
	fL4[1][0] = fGradDetF[1][1]*fGradDetF[0][0];//dux*dvy
	fL4[1][1] = fGradDetF[1][1]*fGradDetF[0][1];//dvx*dvy
	fL4[1][2] = fGradDetF[1][1]*fGradDetF[0][2];//dwx*dvy
	fL4[2][0] = fGradDetF[1][2]*fGradDetF[0][0];//dux*dwy
	fL4[2][1] = fGradDetF[1][2]*fGradDetF[0][1];//dvx*dwy
	fL4[2][2] = fGradDetF[1][2]*fGradDetF[0][2];//dwx*dwy
	
	fL5[0][0] = fGradDetF[1][0]*fGradDetF[1][0];//duy*duy
	fL5[0][1] = fGradDetF[1][0]*fGradDetF[1][1];//dvy*duy
	fL5[0][2] = fGradDetF[1][0]*fGradDetF[1][2];//dwy*duy
	fL5[1][0] = fGradDetF[1][1]*fGradDetF[1][0];//duy*dvy
	fL5[1][1] = fGradDetF[1][1]*fGradDetF[1][1];//dvy*dvy
	fL5[1][2] = fGradDetF[1][1]*fGradDetF[1][2];//dwy*dvy
	fL5[2][0] = fGradDetF[1][2]*fGradDetF[1][0];//duy*dwy
	fL5[2][1] = fGradDetF[1][2]*fGradDetF[1][1];//dvy*dwy
	fL5[2][2] = fGradDetF[1][2]*fGradDetF[1][2];//dwy*dwy
	
	fL6[0][0] = fGradDetF[1][0]*fGradDetF[2][0];//duz*duy
	fL6[0][1] = fGradDetF[1][0]*fGradDetF[2][1];//dvz*duy
	fL6[0][2] = fGradDetF[1][0]*fGradDetF[2][2];//dwz*duy
	fL6[1][0] = fGradDetF[1][1]*fGradDetF[2][0];//duz*dvy
	fL6[1][1] = fGradDetF[1][1]*fGradDetF[2][1];//dvz*dvy
	fL6[1][2] = fGradDetF[1][1]*fGradDetF[2][2];//dwz*dvy
	fL6[2][0] = fGradDetF[1][2]*fGradDetF[2][0];//duz*dwy
	fL6[2][1] = fGradDetF[1][2]*fGradDetF[2][1];//dvz*dwy
	fL6[2][2] = fGradDetF[1][2]*fGradDetF[2][2];//dwz*dwy
	
	fL7[0][0] = fGradDetF[2][0]*fGradDetF[0][0];//dux*duz
	fL7[0][1] = fGradDetF[2][0]*fGradDetF[0][1];//dvx*duz
	fL7[0][2] = fGradDetF[2][0]*fGradDetF[0][2];//dwx*duz
	fL7[1][0] = fGradDetF[2][1]*fGradDetF[0][0];//dux*dvz
	fL7[1][1] = fGradDetF[2][1]*fGradDetF[0][1];//dvx*dvz
	fL7[1][2] = fGradDetF[2][1]*fGradDetF[0][2];//dwx*dvz
	fL7[2][0] = fGradDetF[2][2]*fGradDetF[0][0];//dux*dwz
	fL7[2][1] = fGradDetF[2][2]*fGradDetF[0][1];//dvx*dwz
	fL7[2][2] = fGradDetF[2][2]*fGradDetF[0][2];//dwx*dwz
	
	fL8[0][0] = fGradDetF[2][0]*fGradDetF[1][0];//duy*duz
	fL8[0][1] = fGradDetF[2][0]*fGradDetF[1][1];//dvy*duz
	fL8[0][2] = fGradDetF[2][0]*fGradDetF[1][2];//dwy*duz
	fL8[1][0] = fGradDetF[2][1]*fGradDetF[1][0];//duy*dvz
	fL8[1][1] = fGradDetF[2][1]*fGradDetF[1][1];//dvy*dvz
	fL8[1][2] = fGradDetF[2][1]*fGradDetF[1][2];//dwy*dvz
	fL8[2][0] = fGradDetF[2][2]*fGradDetF[1][0];//duy*dwz
	fL8[2][1] = fGradDetF[2][2]*fGradDetF[1][1];//dvy*dwz
	fL8[2][2] = fGradDetF[2][2]*fGradDetF[1][2];//dwy*dwz
	
	fL9[0][0] = fGradDetF[2][0]*fGradDetF[2][0];//duz*duz
	fL9[0][1] = fGradDetF[2][0]*fGradDetF[2][1];//dvz*duz
	fL9[0][2] = fGradDetF[2][0]*fGradDetF[2][2];//dwz*duz
	fL9[1][0] = fGradDetF[2][1]*fGradDetF[2][0];//duz*dvz
	fL9[1][1] = fGradDetF[2][1]*fGradDetF[2][1];//dvz*dvz
	fL9[1][2] = fGradDetF[2][1]*fGradDetF[2][2];//dwz*dvz
	fL9[2][0] = fGradDetF[2][2]*fGradDetF[2][0];//duz*dwz
	fL9[2][1] = fGradDetF[2][2]*fGradDetF[2][1];//dvz*dwz
	fL9[2][2] = fGradDetF[2][2]*fGradDetF[2][2];//dwz*dwz
	
	
	REAL detF = (ux+1.)*(vy+1.)*(wz+1.) + vx*wy*uz + wx*uy*vz - wx*(vy+1.)*uz - (ux+1.)*wy*vz - vx*uy*(wz+1.);
	
	if(detF < 0) {
		cout << "\nDeterminante negativo!\n";
	}
	
	REAL fac1 = 2.*fCoef1 - fCoef2/detF/detF;//A
	REAL fac2 = 2.*fCoef1*detF + fCoef2/detF;//B
	REAL c = fCoef3;//C
	
	int j;
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			global[i][j][0] = fac1*fL1[i][j];
			global[i][j][1] = fac1*fL2[i][j]+fac2*fK2[i][j];
			global[i][j][2] = fac1*fL3[i][j]+fac2*fK3[i][j];
			global[i][j][3] = fac1*fL4[i][j]+fac2*fK4[i][j];
			global[i][j][4] = fac1*fL5[i][j];
			global[i][j][5] = fac1*fL6[i][j]+fac2*fK6[i][j];
			global[i][j][6] = fac1*fL7[i][j]+fac2*fK7[i][j];
			global[i][j][7] = fac1*fL8[i][j]+fac2*fK8[i][j];
			global[i][j][8] = fac1*fL9[i][j];
		}
		global[i][i][0] += c*fE1[i];
		global[i][i][4] += c*fE5[i];
		global[i][i][8] += c*fE9[i];
	}
	
	//AT�AQUI AS TRES MATRIZES A INTEGRAR S� fK , fL E fE , SENDO fE CONSTANTE
	REAL gradJx[3],gradJy[3],gradJz[3];
	REAL gradtrCx[3],gradtrCy[3],gradtrCz[3];
	int k,l;
	int nshape = phi.Rows();
	for(i=0; i<3; i++) {
		gradJx[i] = fGradDetF[0][i];//dux , dvx , dwx
		gradJy[i] = fGradDetF[1][i];//duy , dvy , dwy
		gradJz[i] = fGradDetF[2][i];//duz , dvz , dwz
		gradtrCx[i] = fGradtrC[0][i];//dux , dvx , dwx
		gradtrCy[i] = fGradtrC[1][i];//duy , dvy , dwy
		gradtrCz[i] = fGradtrC[2][i];//duz , dvz , dwz
	}
	
	REAL *efptr = &ef(0,0);
	REAL *dphiptr = &dphi(0,0);
	int nrowek = ek.Rows();
	REAL *ekptr = &ek(0,0);
	for(k=0; k<3; k++) {
		REAL kval[3] = {(fac2*gradJx[k]+c*gradtrCx[k]),(fac2*gradJy[k]+c*gradtrCy[k]),(fac2*gradJz[k]+c*gradtrCz[k])};
		for(i=0; i<nshape; i++) {
			efptr[i*3+k] += -weight*(dphiptr[3*i]*kval[0]+
									 dphiptr[3*i+1]*kval[1]+
									 dphiptr[3*i+2]*kval[2]);//G = grad W
		}
	}
	for(k=0; k<3; k++) {
		for(l=0; l<3; l++) {
			for(i=0; i<nshape; i++) {
				for(j=0; j<nshape; j++) {
					ekptr[i*3+k+nrowek*(j*3+l)] += weight*(
														   dphiptr[3*i]*dphiptr[3*j]*global[k][l][0] + dphiptr[3*i]*dphiptr[1+3*j]*global[k][l][1] +
														   dphiptr[3*i]*dphiptr[2+3*j]*global[k][l][2] + dphiptr[1+3*i]*dphiptr[3*j]*global[k][l][3] +
														   dphiptr[1+3*i]*dphiptr[1+3*j]*global[k][l][4] + dphiptr[1+3*i]*dphiptr[2+3*j]*global[k][l][5] +
														   dphiptr[2+3*i]*dphiptr[3*j]*global[k][l][6] + dphiptr[2+3*i]*dphiptr[1+3*j]*global[k][l][7] +
														   dphiptr[2+3*i]*dphiptr[2+3*j]*global[k][l][8] );
				}
			}
		}
	}
	
#endif
}
/*for(idf=0; idf<3; idf++) {
 for(jdf=0; jdf<3; jdf++)
 ek(i*3+idf,j*3+jdf) += KIJ(idf,jdf);
 }*/


/** returns the variable index associated with the name*/
int TPZMatHyperElastic::VariableIndex(const std::string &name){
	if(!strcmp("Displacement6",name.c_str()))   return  1;
	if(!strcmp("displacement",name.c_str()))     return  2;
	if(!strcmp("Solution",name.c_str()))     return  2;
	if(!strcmp("Derivate",name.c_str()))     return  3;
	if(!strcmp("VonMises",name.c_str())) return 4;
	if(!strcmp("POrder",name.c_str()))       return 10;
	return TPZMaterial::VariableIndex(name);
}

int TPZMatHyperElastic::NSolutionVariables(int var){
	
	if(var == 1) return 6;
	if(var == 2) return 3;
	if(var == 3) return 9;
	if(var == 4) return 1;
	if(var == 10) return 1;
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatHyperElastic::Solution(TPZVec<REAL> &Sol,TPZFMatrix<REAL> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<REAL> &Solout){
	
	if(var == 1) Solout.Resize(6,0.);
	if(var == 2) Solout.Resize(3,0.);
	if(var == 1|| var == 2) {
      	Solout[0] = Sol[0];//function
		Solout[1] = Sol[1];//function
		Solout[2] = Sol[2];//function
	} else
		if(var == 3) {
			Solout.Resize(9);
			int k=0;
			for(int i=0;i<3;i++) {
				Solout[k++] = DSol(0,i);//derivate
				Solout[k++] = DSol(1,i);//derivate
				Solout[k++] = DSol(2,i);//derivate
			}
		}
		else if(var == 10) {
			Solout.Resize(1);
			Solout[0] = 1.;
		}
		else if(var == 4) {
			
			TPZFMatrix<REAL> F(DSol),Ft(3,3,0.0),I(3,3,0.),S(3,3,0.);
			F(0,0)+=1;
			F(1,1)+=1;
			F(2,2)+=1;//a diagonal e' sempre a mesma
			F.Transpose(&Ft);
			TPZFMatrix<REAL> B = F*Ft;
			
			REAL ux,uy,uz,vx,vy,vz,wx,wy,wz;
			ux = DSol(0,0);
			uy = DSol(1,0);
			uz = DSol(2,0);
			
			vx = DSol(0,1);
			vy = DSol(1,1);
			vz = DSol(2,1);
			
			wx = DSol(0,2);
			wy = DSol(1,2);
			wz = DSol(2,2);
			
			REAL J = (ux+1.)*(vy+1.)*(wz+1.) + vx*wy*uz + wx*uy*vz - wx*(vy+1.)*uz - (ux+1.)*wy*vz - vx*uy*(wz+1.);
			
			I(0,0)=1.0;
			I(1,1)=1.0;
			I(2,2)=1.0;
			TPZFMatrix<REAL> sigmaF = (REAL(fNu/J)*B+REAL((fLambda*REAL(0.5)*(J*J-1.0)-fNu)/J)*I);
			REAL trsigma = sigmaF(0,0)+ sigmaF(1,1)+ sigmaF(2,2);
			S = sigmaF - REAL(trsigma/3.0)*I;
			int i,j;
			REAL J2=0.0;
			for(i=0;i<3;i++)  for(j=0;j<3;j++) J2 += S(i,j)* S(i,j);
			Solout[0] = sqrt(3.0*J2);
			
			
		}
		else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TPZMatHyperElastic::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix<REAL> &/*DSol*/, TPZFMatrix<REAL> &/*axes*/, TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix<REAL> &DSol, TPZFMatrix<REAL> &axes, TPZVec<REAL> &flux)
}

void TPZMatHyperElastic::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
								TPZFMatrix<REAL> &dudx, TPZFMatrix<REAL> &axes, TPZVec<REAL> &/*flux*/,
								TPZVec<REAL> &u_exact,TPZFMatrix<REAL> &du_exact,TPZVec<REAL> &values) {
	
	//TPZVec<REAL> sol(1),dsol(3);
	TPZManVector<REAL> sol(3),dsol(9);
	Solution(u,dudx,axes,2,sol);
	Solution(u,dudx,axes,3,dsol);
	//values[1] : erro em norma L2
	values[1]  = pow(sol[0] - u_exact[0],(REAL)2.0);
	values[1] += pow(sol[1] - u_exact[1],(REAL)2.0);
	values[1] += pow(sol[2] - u_exact[2],(REAL)2.0);
	//values[2] : erro em semi norma H1
	int k=0;
	values[2] = 0.;
	for(int i=0;i<3;i++) {
		values[2] += pow(dsol[k++] - du_exact(0,i),(REAL)2.0);
		values[2] += pow(dsol[k++] - du_exact(1,i),(REAL)2.0);
		values[2] += pow(dsol[k++] - du_exact(2,i),(REAL)2.0);
	}
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}



void TPZMatHyperElastic::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<REAL> &ek, TPZFMatrix<REAL> &ef, TPZBndCond &bc){
	
	//   TPZFMatrix<REAL> &dphi = data.dphix;
	// TPZFMatrix<REAL> &dphiL = data.dphixl;
	// TPZFMatrix<REAL> &dphiR = data.dphixr;
	TPZFMatrix<REAL> &phi = data.phi;
	// TPZFMatrix<REAL> &phiL = data.phil;
	// TPZFMatrix<REAL> &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZVec<REAL> &sol=data.sol[0];
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	//TPZFMatrix<REAL> &dsol=data.dsol;
	// TPZFMatrix<REAL> &dsolL=data.dsoll;
	// TPZFMatrix<REAL> &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix<REAL> &daxesdksi=data.daxesdksi;
	// TPZFMatrix<REAL> &axes=data.axes;
	
	if(bc.Material().operator ->() != this){
		PZError << "TPZMatHyperElastic.ContributeBC : this material don't exists \n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "ContributeBC.aplybc, unknown boundary condition type : "<<bc.Type() << endl;
	}
	
	int ndof = NStateVariables();
	int nnod = ek.Rows()/ndof;
	int r = ndof;
	
	int idf,jdf,in,jn;
	switch(bc.Type()){
		case 0:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += gBigNumber*phi(in,0)*(bc.Val2()(idf,0)-sol[idf])*weight;
				}
				for(jn=0 ; jn<nnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
					}
				}
			}
			break;
			
		case 1:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					//(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0)-sol[idf]);
					(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0));
				}
			}
			break;
			
		case 2:
			for(in=0 ; in<nnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					for (jdf=0; jdf<r; jdf++){
						(ef)(in*r+idf,0) += phi(in,0)*bc.Val1()(idf,jdf)*(bc.Val2()(jdf,0)-sol[jdf])*weight;
					}
					for(jn=0 ; jn<nnod ; ++jn) {
						for(idf = 0;idf<r;idf++) {
							for(jdf = 0;jdf<r;jdf++) {
								ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
							}
						}
					}
				}
				
            }
			
	}//fim switch
}


#ifdef _AUTODIFF


/* The function below makes the correspondence between the dsol vector
 and a matrix ordered F operator
 */
inline int ith(const int i, const int j)
{
	return i*3+j;
}


void TPZMatHyperElastic::ContributeEnergy(TPZVec<REAL> &x,
										  TPZVec<FADFADREAL> &sol, TPZVec<FADFADREAL> &dsol,
										  FADFADREAL &U, REAL weight)
{
	//     const int dim = 3;
	//  cout <<  "Solution gradient\n ";
	//  int i;
	//  for(i=0; i<9; i++) cout  << dsol[i].val().val() <<  " ";
	//  cout << endl;
	FADFADREAL J, TrC; // J = det(F); TrC = Trace(C)
	/*     int numder = dsol[0].size();
     FADFADREAL F[3][3], J, J2, TrC; // J = det(F); TrC = Trace(C)
     FADREAL defaultFAD(numder, 0., 0.);
	 
	 
     FADFADREAL ADOne(numder, FADREAL(numder, 1., 0.), defaultFAD);
     FADFADREAL ADTwo(numder, FADREAL(numder, 2., 0.), defaultFAD);
     FADFADREAL ADThree(numder, FADREAL(numder, 3., 0.), defaultFAD);
     FADFADREAL ADFour(numder, FADREAL(numder, 4., 0.), defaultFAD);
     FADFADREAL ADLambda(numder, FADREAL(numder, fLambda, 0.), defaultFAD);
     FADFADREAL ADNu(numder, FADREAL(numder, fNu, 0.), defaultFAD);
     TrC = FADFADREAL(numder, defaultFAD, defaultFAD);
	 
     int i,j,k;
     for(i = 0; i < dim; i++)
     {
	 for(j = 0; j < dim; j++)
	 F[i][j] = dsol[i * dim + j]; // Fij = dUi/dXj // GradU
	 F[i][i] += ADOne; // F = GradU + I
     }
	 
	 // C = Transpose[F] * F
	 // evaluating the diagonal of the Cauchy tensor
	 
	 
     for(i = 0; i < dim; i++)
     {
	 for(k = 0; k < dim; k++)
	 TrC += F[k][i] * F[k][i]; // matrix multiplication
     }
	 
     J = F[0][0] * F[1][1] * F[2][2] +
	 F[0][1] * F[1][2] * F[2][0] +
	 F[0][2] * F[1][0] * F[2][1] -
	 F[0][2] * F[1][1] * F[2][0] -
	 F[0][1] * F[1][0] * F[2][2] -
	 F[0][0] * F[1][2] * F[2][1]; //  J = det(F)
	 
     TPZVecFADFADREAL> DiagF(3);
     DiagF[0] = dsol[    0] + FADREAL(1.);// element [0][0]
     DiagF[1] = dsol[  3+1] + FADREAL(1.);// element [1][1]
     DiagF[2] = dsol[2*3+2] + FADREAL(1.);// element [2][2]
	 
     TrC = DiagF[0]*DiagF[0]+
	 dsol[ith(1,0)] * dsol[ith(1,0)] +
	 dsol[ith(2,0)] * dsol[ith(2,0)] +
	 
	 dsol[ith(0,1)] * dsol[ith(0,1)] +
	 DiagF[1]*DiagF[1]+
	 dsol[ith(2,1)] * dsol[ith(2,1)] +
	 
	 dsol[ith(0,2)] * dsol[ith(0,2)] +
	 dsol[ith(1,2)] * dsol[ith(1,2)] +
	 DiagF[2]*DiagF[2];
	 */
	FADFADREAL DiagF0(dsol[    0]);
	DiagF0.val().val() += 1.;// element [0][0]
	
	FADFADREAL DiagF1(dsol[  3+1]);
	DiagF1.val().val() += 1.;// element [0][0]
	
	FADFADREAL DiagF2(dsol[2*3+2]);
	DiagF2.val().val() += 1.;// element [0][0]
	/*
     TrC = DiagF0*DiagF0 +
	 dsol[ith(1,0)] * dsol[ith(1,0)] +
	 dsol[ith(2,0)] * dsol[ith(2,0)] +
	 
	 dsol[ith(0,1)] * dsol[ith(0,1)] +
	 DiagF1*DiagF1 +
	 dsol[ith(2,1)] * dsol[ith(2,1)] +
	 
	 dsol[ith(0,2)] * dsol[ith(0,2)] +
	 dsol[ith(1,2)] * dsol[ith(1,2)] +
	 DiagF2*DiagF2;
	 
     J = DiagF0          * DiagF1         * DiagF2 +
	 dsol[ith(0,1)] * dsol[ith(1,2)] * dsol[ith(2,0)] +
	 dsol[ith(0,2)] * dsol[ith(1,0)] * dsol[ith(2,1)] -
	 dsol[ith(0,2)] * DiagF1         * dsol[ith(2,0)] -
	 dsol[ith(0,1)] * dsol[ith(1,0)] * DiagF2 -
	 DiagF0         * dsol[ith(1,2)] * dsol[ith(2,1)]; //  J = det(F)
	 
     U += (J*J - FADREAL(1.)) * FADREAL(fLambda/4.) -
	 log( J ) * FADREAL(fLambda/2.+fNu) +
	 (TrC - FADREAL(3.)) * FADREAL(fNu/2.);
	 
	 */
	TrC =  DiagF0*DiagF0;
	TrC += dsol[ith(1,0)] * dsol[ith(1,0)];
	TrC += dsol[ith(2,0)] * dsol[ith(2,0)];
	
	TrC += dsol[ith(0,1)] * dsol[ith(0,1)];
	TrC += DiagF1*DiagF1;
	TrC += dsol[ith(2,1)] * dsol[ith(2,1)];
	
	TrC += dsol[ith(0,2)] * dsol[ith(0,2)];
	TrC += dsol[ith(1,2)] * dsol[ith(1,2)];
	TrC += DiagF2*DiagF2;
	
	//     cout <<  "TrC\n" << TrC << endl;
	J = DiagF0          * DiagF1         * DiagF2;
	J += dsol[ith(0,1)] * dsol[ith(1,2)] * dsol[ith(2,0)];
	J += dsol[ith(0,2)] * dsol[ith(1,0)] * dsol[ith(2,1)];
	J -= dsol[ith(0,2)] * DiagF1         * dsol[ith(2,0)];
	J -= dsol[ith(0,1)] * dsol[ith(1,0)] * DiagF2;
	J -= DiagF0         * dsol[ith(1,2)] * dsol[ith(2,1)]; //  J = det(F)
	
	//     cout <<  "J\n " << J << endl;
	U += (J*J - FADREAL(1.)) * FADREAL(weight*fLambda/4.);
	U -= log( J ) * FADREAL(weight*(fLambda/2.+fNu));
	U += (TrC - FADREAL(3.)) * FADREAL(weight*fNu/2.);
	
	//     cout << U << endl;
}

void TPZMatHyperElastic::ContributeBCEnergy(TPZVec<REAL> & x,
											TPZVec<FADFADREAL> & sol, FADFADREAL &U,
											REAL weight, TPZBndCond &bc)
{
	if(bc.Material().operator->() != this){
		PZError << "TPZMatHyperElastic.ContributeBC : this material doesn't exist \n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "ContributeBC.aplybc, unknown boundary condition type : "<<bc.Type() << endl;
	}
	
	TPZVec<FADFADREAL> solMinBC(3), BCsolMinBC(3);
	solMinBC[0] = sol[0] - FADREAL(bc.Val2()(0,0));
	solMinBC[1] = sol[1] - FADREAL(bc.Val2()(1,0));
	solMinBC[2] = sol[2] - FADREAL(bc.Val2()(2,0));
	
	switch(bc.Type()){
		case 0:// Dirichlet condition
            // U += 1/2* Big * weight * Integral((u - u0)^2 dOmega)
			U += ( (solMinBC[0] * solMinBC[0]) +
				  (solMinBC[1] * solMinBC[1]) +
				  (solMinBC[2] * solMinBC[2]) ) *
			FADREAL(gBigNumber * weight / 2.);
			break;
		case 1:// Neumann condition
            // U -= weight * Integral([g].u dOmega)
			U -= ( sol[0] * FADREAL(bc.Val2()(0,0) ) +
				  sol[1] * FADREAL(bc.Val2()(1,0) ) +
				  sol[2] * FADREAL(bc.Val2()(2,0) ) ) *
			FADREAL(weight);
			break;
		case 2:// condi�o mista
            // U += 1/2 * weight * Integral(<(u-u0), [g].(u-u0)> dOmega)
			BCsolMinBC[0] = solMinBC[0] * FADREAL(bc.Val1()(0,0)) +
			solMinBC[1] * FADREAL(bc.Val1()(0,1)) +
			solMinBC[2] * FADREAL(bc.Val1()(0,2));
			BCsolMinBC[1] = solMinBC[0] * FADREAL(bc.Val1()(1,0)) +
			solMinBC[1] * FADREAL(bc.Val1()(1,1)) +
			solMinBC[2] * FADREAL(bc.Val1()(1,2));
			BCsolMinBC[2] = solMinBC[0] * FADREAL(bc.Val1()(2,0)) +
			solMinBC[1] * FADREAL(bc.Val1()(2,1)) +
			solMinBC[2] * FADREAL(bc.Val1()(2,2));
			U += ( solMinBC[0] * BCsolMinBC[0] +
				  solMinBC[1] * BCsolMinBC[1] +
				  solMinBC[2] * BCsolMinBC[2] ) *
			FADREAL(weight / 2.);
			break;
	}
}

void TPZMatHyperElastic::ComputeEnergy(REAL lambda, REAL mu,  TPZFMatrix<REAL> &dsol, TFad<9,TFad<9,REAL> > &energy){
	TFad<9,TFad<9,REAL> > tensor[3][3],J,TrC;
	int i,j;
	//  dsol.Print( "Solution gradient ");
	for(i=0; i<3; i++) {
		for(j=0; j<3; j++) {
			tensor[i][j].val().val() = dsol(j,i);
			tensor[i][j].fastAccessDx(j*3+i).val() = 1.;
			tensor[i][j].val().fastAccessDx(j*3+i) = 1.;
		}
		tensor[i][i].val().val() += 1.;
	}
	TrC = tensor[0][0]*tensor[0][0]+tensor[0][1]*tensor[0][1]+tensor[0][2]*tensor[0][2]+tensor[1][0]*tensor[1][0]+tensor[1][1]*tensor[1][1]+tensor[1][2]*tensor[1][2]+
    tensor[2][0]*tensor[2][0]+tensor[2][1]*tensor[2][1]+tensor[2][2]*tensor[2][2];
	
	//  cout <<  "TrC\n " << TrC << endl;
	J = tensor[0][0] * tensor[1][1] * tensor[2][2] + 
    tensor[0][1] * tensor[1][2] * tensor[2][0] +
    tensor[0][2] * tensor[1][0] * tensor[2][1] -
    tensor[0][2] * tensor[1][1] * tensor[2][0] -
    tensor[0][1] * tensor[1][0] * tensor[2][2] -
    tensor[0][0] * tensor[1][2] * tensor[2][1]; //  J = det(F)
	//  cout <<  "J\n " << J << endl;
	energy = (J*J - TFad<9,REAL>(1.)) * TFad<9,REAL>(lambda/4.) -
    log( J ) * TFad<9,REAL>((lambda/2.+mu)) +
    (TrC - TFad<9,REAL>(3.)) * TFad<9,REAL>(mu/2.);
	
}


#endif
