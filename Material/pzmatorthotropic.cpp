/**
 * \file
 * @brief Contains implementations of the TPZMatOrthotropic methods.
 */

#include "pzmatorthotropic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
//#include "pztempmat.h"
#include "pzmanvector.h"
#include <math.h>
#include <cmath>
#include <fstream>
using namespace std;

TPZMatOrthotropic::TPZMatOrthotropic(int nummat,TPZFMatrix naxes,REAL eppx,REAL eppy,
                                     REAL eppz,REAL vxy,REAL vyz,REAL vzx,
									 REAL gxy,REAL gyz,REAL gzx) :
TPZMaterial(nummat),
fKXX(3,3,0.),fKYY(3,3,0.),fKZZ(3,3,0.),
fKXY(3,3,0.),fKYX(3,3,0.),fKXZ(3,3,0.),
fKZX(3,3,0.),fKYZ(3,3,0.),
fKZY(3,3,0.),fXf(3,1,0.)
{
	/**eixos locais por linhas*/
	//normaliza os eixos
	Normalize(naxes);
	fLocAxs = naxes;
	fGxy = gxy;
	fGyz = gyz;
	fGzx = gzx;
	fEppx = eppx;
	fEppy = eppy;
	fEppz = eppz;
	fVxy  = vxy;
	fVyx  = fEppy*fVxy/fEppx;
	fVyz  = vyz;
	fVzy  = fEppz*fVyz/fEppy;
	fVzx  = vzx;
	fVxz  = fEppx*fVzx/fEppz;
	fNumNom = -1.0 + fVxy*fVyx + fVxz*fVzx  + fVyz*fVzy + fVxy*fVyz*fVzx + fVyx*fVzy*fVxz;
	
	int i,j;
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKXX(i,j) = 0.;   
	fKXX(0,0) = -(fEppx - fEppx*fVyz*fVzy)/fNumNom;
	fKXX(1,1) = fGxy;
	fKXX(2,2) = fGzx;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKXY(i,j) = 0.;//XY 
	fKXY(0,1) = -(fEppx*fVyx + fEppx*fVyz*fVzx)/fNumNom;
	fKXY(1,0) = fGxy;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKXZ(i,j) = 0.;//XZ
	fKXZ(0,2) = -(fEppx*fVzx + fEppx*fVyx*fVzy)/fNumNom;
	fKXZ(2,0) = fGzx;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKYX(i,j) = 0.;//YX
	fKYX(1,0) = -(fEppy*fVxy + fEppy*fVxz*fVzy)/fNumNom;
	fKYX(0,1) = fGxy;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKYY(i,j) = 0.;
	fKYY(0,0) = fGxy; 
	fKYY(1,1) = -(fEppy - fEppy*fVxz*fVzx)/fNumNom;
	fKYY(2,2) = fGyz;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKYZ(i,j) = 0.;//YZ
	fKYZ(1,2) = -(fEppy*fVzy + fEppy*fVxy*fVzx)/fNumNom;
	fKYZ(2,1) = fGyz;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKZX(i,j) = 0.;//ZX
	fKZX(2,0) = -fEppz*(fVxz + fVxy*fVyz)/fNumNom;
	fKZX(0,2) = fGzx;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKZY(i,j) = 0.;//ZY
	fKZY(2,1) = -fEppz*(fVyz + fVxz*fVyx)/fNumNom;
	fKZY(1,2) = fGyz;
	
	for(i=0; i<3; i++) for(j=0; j<3; j++) fKZZ(i,j) = 0.;   
	fKZZ(0,0) = fGzx;
	fKZZ(1,1) = fGyz;
	fKZZ(2,2) = -(fEppz - fEppz*fVxy*fVyx)/fNumNom;
}

TPZMatOrthotropic::~TPZMatOrthotropic() {
}

int TPZMatOrthotropic::NStateVariables() {
	return 3;
}

void TPZMatOrthotropic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "Eppx = " << fEppx << "\tEppy = " << fEppy << "\tEppz = " << fEppz << endl;
	out << "vxy = " << fVxy << "\tvyx = " << fVyx << endl;
	out << "vzx = " << fVzx << "\tvxz = " << fVxz << endl;
	out << "vxy = " << fVxy << "\tvyx = " << fVyx << endl;
	out << "gxy = " << fGxy << "\tgyz = " << fGyz << "\tgzx = " << fGzx << endl;
	out << "numnomvar = " << fNumNom << endl; 
	out << "\n>>>>>>>>>>>>>> MATRIZES KIJ <<<<<<<<<<<<<<<<<<\n\n";
	fKXX.Print("KXX",out);  
	fKXY.Print("KXY",out);  
	fKXZ.Print("KXZ",out);  
	fKYX.Print("KYX",out);  
	fKYY.Print("KYY",out);  
	fKYZ.Print("KYZ",out);  
	fKZX.Print("KZX",out);  
	fKZY.Print("KZY",out);  
	fKZZ.Print("KZZ",out);  
	
	TPZMaterial::Print(out);
}

/** @brief Output file to write stiffness matrix */
ofstream MatrizesK("MatrizesK.out");
void TPZMatOrthotropic::Contribute(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix &ek,
                                   TPZFMatrix &ef) {
	
	TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	TPZFMatrix &axes=data.axes;
	
	int phr = phi.Rows();
	if(fForcingFunction) {            
		TPZManVector<REAL> res(3);//,&fXf(0,0),3);// phi(in, 0) = phi_in
		fForcingFunction->Execute(x,res);    // dphi(i,j) = dphi_j/dxi
		int i;
		for(i=0; i<3; i++) fXf(i,0) = res[i];
	}
	
	TPZFMatrix Rt(fLocAxs),R(3,3,0.),newaxes(3,3,0.);
	TPZFMatrix kxx(3,3,0.),kxy(3,3,0.),kxz(3,3,0.);
	TPZFMatrix kyx(3,3,0.),kyy(3,3,0.),kyz(3,3,0.);
	TPZFMatrix kzx(3,3,0.),kzy(3,3,0.),kzz(3,3,0.);
	
	Rt.Transpose(&R);
	/** as linhas de axes s� as derivadas nas direc�s do jacobiano,
     * cada linha de newaxes �o produto interno entre uma direc� (ni)
     * da placa e as 3 direc�s locais da derivada (a0,a1,a2), assim a
     * linha i de newaxes �  ni*a0 ni*a1 ni*a2
     */
	axes.Transpose(&newaxes);
	newaxes = Rt*newaxes;
	
	kxx = R*(fKXX*Rt); 
	kxy = R*(fKXY*Rt); 
	kxz = R*(fKXZ*Rt);
	kyx = R*(fKYX*Rt); 
	kyy = R*(fKYY*Rt); 
	kyz = R*(fKYZ*Rt); 
	kzx = R*(fKZX*Rt); 
	kzy = R*(fKZY*Rt); 
	kzz = R*(fKZZ*Rt); 
	
	MatrizesK << "\n>>>>>>>>>>>>>> MATRIZES KIJ <<<<<<<<<<<<<<<<<<\n\n";
	kxx.Print("KXX",MatrizesK);  
	kxy.Print("KXY",MatrizesK);  
	kxz.Print("KXZ",MatrizesK);  
	kyx.Print("KYX",MatrizesK);  
	kyy.Print("KYY",MatrizesK);  
	kyz.Print("KYZ",MatrizesK);  
	kzx.Print("KZX",MatrizesK);  
	kzy.Print("KZY",MatrizesK);  
	kzz.Print("KZZ",MatrizesK);  
	
	for( int in = 0; in < phr; in++ ) {
		ef(3*in, 0)   += weight * phi(in, 0) * fXf(0,0);
		ef(3*in+1, 0) += weight * phi(in, 0) * fXf(1,0);
		ef(3*in+2, 0) += weight * phi(in, 0) * fXf(2,0);
		for( int jn = 0; jn < phr; jn++ ) {
			
			/**derivadas com respeito a in*/         
			REAL dphin0 = newaxes(0,0)*dphi(0,in)+newaxes(0,1)*dphi(1,in)+newaxes(0,2)*dphi(2,in);/**produto a0*n0, a1*n0, a2*n0*/
			REAL dphin1 = newaxes(1,0)*dphi(0,in)+newaxes(1,1)*dphi(1,in)+newaxes(1,2)*dphi(2,in);/**produto a0*n1, a1*n1, a2*n1*/
			REAL dphin2 = newaxes(2,0)*dphi(0,in)+newaxes(2,1)*dphi(1,in)+newaxes(2,2)*dphi(2,in);/**produto a0*n2, a1*n2, a2*n2*/
			/**derivadas com respeito a jn*/
			REAL dphjn0 = newaxes(0,0)*dphi(0,jn)+newaxes(0,1)*dphi(1,jn)+newaxes(0,2)*dphi(2,jn);/**produto a0*n0, a1*n0, a2*n0*/
			REAL dphjn1 = newaxes(1,0)*dphi(0,jn)+newaxes(1,1)*dphi(1,jn)+newaxes(1,2)*dphi(2,jn);/**produto a0*n1, a1*n1, a2*n1*/
			REAL dphjn2 = newaxes(2,0)*dphi(0,jn)+newaxes(2,1)*dphi(1,jn)+newaxes(2,2)*dphi(2,jn);/**produto a0*n2, a1*n2, a2*n2*/
			
			int k,l;
			for(k=0; k<3; k++) {
				for(l=0; l<3; l++) {
					ek(3*in+k,3*jn+l) += weight * ( dphin0*dphjn0*kxx(k,l)  + 
												   dphin0*dphjn1*kxy(k,l)  + 
												   dphin0*dphjn2*kxz(k,l)  +
												   dphin1*dphjn0*kyx(k,l)  + 
												   dphin1*dphjn1*kyy(k,l)  + 
												   dphin1*dphjn2*kyz(k,l)  +
												   dphin2*dphjn0*kzx(k,l)  + 
												   dphin2*dphjn1*kzy(k,l)  + 
												   dphin2*dphjn2*kzz(k,l)
												   );
				}
			}
		}
	}
}

void TPZMatOrthotropic::ContributeBC(TPZMaterialData &data,
                                     REAL weight,
                                     TPZFMatrix &ek,
                                     TPZFMatrix &ef,
                                     TPZBndCond &bc) {
	
	// TPZFMatrix &dphi = data.dphix;
	// TPZFMatrix &dphiL = data.dphixl;
	// TPZFMatrix &dphiR = data.dphixr;
	TPZFMatrix &phi = data.phi;
	// TPZFMatrix &phiL = data.phil;
	// TPZFMatrix &phiR = data.phir;
	// TPZManVector<REAL,3> &normal = data.normal;
	// TPZManVector<REAL,3> &x = data.x;
	// int &POrder=data.p;
	// int &LeftPOrder=data.leftp;
	// int &RightPOrder=data.rightp;
	// TPZVec<REAL> &sol=data.sol;
	// TPZVec<REAL> &solL=data.soll;
	// TPZVec<REAL> &solR=data.solr;
	// TPZFMatrix &dsol=data.dsol;
	// TPZFMatrix &dsolL=data.dsoll;
	// TPZFMatrix &dsolR=data.dsolr;
	// REAL &faceSize=data.HSize;
	// TPZFMatrix &daxesdksi=data.daxesdksi;
	// TPZFMatrix &axes=data.axes;
	
	const REAL BIGNUMBER  = 1.e12;
	
	int phr = phi.Rows();
	short in,jn,k,l;
	REAL v2[3];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	v2[2] = bc.Val2()(2,0);
	
	switch (bc.Type()) {
		case 0 :// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(3*in,0)   += BIGNUMBER * v2[0] * phi(in,0) * weight;// forced v2[0] x displacement
				ef(3*in+1,0) += BIGNUMBER * v2[1] * phi(in,0) * weight;// forced v2[1] y displacement
				ef(3*in+2,0) += BIGNUMBER * v2[2] * phi(in,0) * weight;// forced v2[2] y displacement
				for (jn = 0 ; jn < phi.Rows(); jn++) {		
					for(k=0; k<3; k++) {
						ek(3*in+k,3*jn+k) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
					}//o Big. deve ser s�na diagonal: OK!
				}
			}          
			break;
			
		case 1 :// Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {// componentes da tra�o normal ao contorno
				ef(3*in,0)   += v2[0] * phi(in,0) * weight; // tra�o em x (ou press�)
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight; // tra�o em y (ou press�)
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight; // tra�o em z (ou press�)
			}
			break;
			
		case 2 :// condi�o mista
			for(in = 0 ; in < phi.Rows(); in++) {          //Sigmaij
				ef(3*in, 0)   += v2[0] * phi(in, 0) * weight;// Neumann x
				ef(3*in+1, 0) += v2[1] * phi(in, 0) * weight;// Neumann y
				ef(3*in+2, 0) += v2[2] * phi(in, 0) * weight;// Neumann z
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					for(k=0; k<3; k++) {
						for(l=0; l<3; l++) {
							ek(3*in+k,3*jn+l) += bc.Val1()(k,l) * phi(in,0) * phi(jn,0) * weight;//integral de contorno
						}
					}
				}
			}
	}
}


/** returns the variable index associated with the name*/
int TPZMatOrthotropic::VariableIndex(const std::string &name){
	if("DisplacX" == name) return 0;
	if(!strcmp("DisplacX",name.c_str()))          return  0;
	if(!strcmp("DisplacY",name.c_str()))          return  1;
	if(!strcmp("DisplacZ",name.c_str()))          return  2;
	if(!strcmp("DisplacN0",name.c_str()))         return  3;
	if(!strcmp("DisplacN1",name.c_str()))         return  4;
	if(!strcmp("DisplacN2",name.c_str()))         return  5;
	if(!strcmp("GlobDisplac",name.c_str()))       return  6;
	if(!strcmp("FiberDisplac",name.c_str()))      return  7;
	if(!strcmp("Tension",name.c_str()))           return  8;
	if(!strcmp("Tensor",name.c_str()))            return  9;
	if(!strcmp("SigX",name.c_str()))              return 10;
	if(!strcmp("SigY",name.c_str()))              return 11;
	if(!strcmp("SigZ",name.c_str()))              return 12;
	if(!strcmp("TauXY",name.c_str()))             return 13;
	if(!strcmp("TauXZ",name.c_str()))             return 14;
	if(!strcmp("TauYZ",name.c_str()))             return 15;
	
	
	cout << "TPZMatOrthotropic::VariableIndex Error\n";
	return -1;
}

int TPZMatOrthotropic::NSolutionVariables(int var){
	
	if(var > -1 && var < 6)  return 1;//escalares
	if(var == 6 || var == 7) return 3;//vetores
	if(var == 8)             return 6;//deve ser
	if(var == 9)             return 9;//tensor completo (sim�rico)
	if(var > 9 && var < 16)  return 1;
	cout << "TPZMatOrthotropic::NSolutionVariables Error\n";
	return 0;
}

void TPZMatOrthotropic::Solution(TPZVec<REAL> &Sol,TPZFMatrix &DSol,TPZFMatrix &axes,int var,TPZVec<REAL> &Solout){ //OBS.:acrescentado ostream
	
	if(var == 0){
		Solout.Resize(1);
		Solout[0] = Sol[0];//eixo global X
		return;
	}
	if(var == 1){
		Solout.Resize(1);
		Solout[0] = Sol[1];//eixo global Y
		return;
	}
	if(var == 2){
		Solout.Resize(1);
		Solout[0] = Sol[2];//eixo global Z
		return;
	}
	if(var == 6){
		Solout.Resize(3);
		Solout[0] = Sol[0];//eixos globais
		Solout[1] = Sol[1];// { X , Y , Z }
		Solout[2] = Sol[2];//
		return;
	}
	if(var==3 || var==4 || var==5 || var==7){  
		cout << "TPZMatOrthotropic::Solution not implemented\n";
		return;
		TPZFMatrix SolN(3,1,0.),sol(3,1,0.);
		sol(0,0) = Sol[0];
		sol(1,0) = Sol[1];
		sol(2,0) = Sol[2];
		SolN = fLocAxs*sol;
		if(var == 3){
			Solout.Resize(1);
			Solout[0] = SolN(0,0);//eixo local n0
			return;
		}
		if(var == 4){
			Solout.Resize(1);
			Solout[0] = SolN(1,0);//eixo local n1
			return;
		}
		if(var == 5){
			Solout.Resize(1);
			Solout[0] = SolN(2,0);//eixo local n2
			return;
		}
		if(var == 7){
			Solout.Resize(3);
			Solout[0] = SolN(0,0);//eixos locais
			Solout[1] = SolN(1,0);// { n0 , n1 , n2 }
			Solout[2] = SolN(2,0);//
			return;
		}
	}
	if(var > 7 && var < 16) {    
		Solout.Resize(6);
		TPZFMatrix axest(3,3,0.),floc_axt(3,3,0.),grdUGdn(3,3,0.),grdUlocdn(3,3,0.);
		axes.Transpose(&axest);
		floc_axt = fLocAxs*axest;
		grdUGdn = floc_axt*DSol;
		TPZFMatrix grdUGdnT(3,3,0.),grdULdn(3,3,0.),grdULdnT(3,3,0.);
		grdUGdn.Transpose(&grdUGdnT);
		grdULdn = fLocAxs*grdUGdnT;
		grdULdn.Transpose(&grdULdnT);
		TPZFMatrix FibStrain(3,3,0.);
		FibStrain = 0.5*(grdULdnT + grdULdn);
		
		REAL epsx  = FibStrain(0,0);// du/dx
		REAL epsy  = FibStrain(1,1);// dv/dy
		REAL epsz  = FibStrain(2,2);// dw/dz
		REAL epsxy = FibStrain(0,1);
		REAL epsyz = FibStrain(1,2);
		REAL epszx = FibStrain(2,0);
		
		REAL SigX = -fEppx*(epsx*(1.-fVyz*fVzy)+epsy*(fVyx+fVyz*fVzx)+epsz*(fVzx+fVyx*fVzy))/fNumNom;
		REAL SigY = -fEppy*(epsy*(1.-fVxz*fVzx)+epsx*(fVxy+fVxz*fVzy)+epsz*(fVzy+fVxy*fVzx))/fNumNom;
		REAL SigZ = -fEppz*(epsz*(1.-fVxy*fVyx)+epsx*(fVxz+fVxy*fVyz)+epsy*(fVyz+fVxz*fVyx))/fNumNom;
		
		REAL TauXY = 2.*fGxy*epsxy;
		REAL TauYZ = 2.*fGyz*epsyz;
		REAL TauZX = 2.*fGzx*epszx;
		
		TPZFMatrix Tensor(3,3,0.);
		Tensor(0,0) = SigX;
		Tensor(1,1) = SigY;
		Tensor(2,2) = SigZ;
		Tensor(0,1) = Tensor(1,0) = TauXY;
		Tensor(1,2) = Tensor(2,1) = TauYZ;
		Tensor(2,0) = Tensor(0,2) = TauZX;
		
		//trecho adicionado apenas para conferencia
		//     cout << "Tensor(0,0) = SigX" << "  =  " << SigX << endl;
		//     cout << "Tensor(1,1) = SigY" << "  =  " << SigY << endl;
		//     cout << "Tensor(2,2) = SigZ" << "  =  " << SigZ << endl;
		//     cout << "Tensor(0,1) = Tensor(1,0) = TauXY" << "  =  " << TauXY << endl;
		//     cout << "Tensor(1,2) = Tensor(2,1) = TauYZ" << "  =  " << TauYZ << endl;
		//     cout << "Tensor(2,0) = Tensor(0,2) = TauZX" << "  =  " << TauZX << endl;
		//final do trecho
		
		TPZFMatrix locaxsT(3,3,0.);
		fLocAxs.Transpose(&locaxsT);
		TPZFMatrix SigGlob = locaxsT*(Tensor*fLocAxs);
		
		if(var == 10){
			Solout.Resize(1);
			Solout[0] = SigX;
			return;
		}
		if(var == 11){
			Solout.Resize(1);
			Solout[0] = SigY;
			return;
		}
		if(var == 12){
			Solout.Resize(1);
			Solout[0] = SigZ;
			return;
		}
		if(var == 13){
			Solout.Resize(1);
			Solout[0] = TauXY;
			return;
		}
		if(var == 14){
			Solout.Resize(1);
			Solout[0] = TauZX;
			return;
		}
		if(var == 15){
			Solout.Resize(1);
			Solout[0] = TauYZ;
			return;
		}
		
		/**ORIGINAL*/
		Solout[0] = SigGlob(0,0);//SigX
		Solout[1] = SigGlob(1,1);//SigY
		Solout[2] = SigGlob(2,2);//SigZ
		Solout[3] = SigGlob(0,1);//TauXY
		Solout[4] = SigGlob(1,2);//TauYZ
		Solout[5] = SigGlob(2,0);//TauZX
		
		if(var == 8) return;
		
		//aqui var �9 e o tensor �sim�rico
		Solout.Resize(9);
		/**para saida no programa do Cafu (DX)*/
		Solout[0] = SigGlob(0,0);//00
		Solout[1] = SigGlob(0,1);//01=10
		Solout[2] = SigGlob(2,0);//02=20
		Solout[3] = SigGlob(0,1);//10=01
		Solout[4] = SigGlob(1,1);//11
		Solout[5] = SigGlob(1,2);//12=21
		Solout[6] = SigGlob(2,0);//20=02
		Solout[7] = SigGlob(1,2);//21=12
		Solout[8] = SigGlob(2,2);//22
		return;
	}
}

void TPZMatOrthotropic::Flux(TPZVec<REAL> &/*x*/, TPZVec<REAL> &/*Sol*/, TPZFMatrix &/*DSol*/, TPZFMatrix &/*axes*/, TPZVec<REAL> &/*flux*/) {
	//Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
}

void TPZMatOrthotropic::Errors(TPZVec<REAL> &/*x*/,TPZVec<REAL> &u,
							   TPZFMatrix &dudx, TPZFMatrix &axes, TPZVec<REAL> &/*flux*/,
							   TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values) {
	
	//TPZVec<REAL> sol(1),dsol(3);
	TPZManVector<REAL> sol(1),dsol(3);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	if(dudx.Rows()<3) {
		REAL dx = du_exact(0,0)*axes(0,0)+du_exact(1,0)*axes(0,1);
		REAL dy = du_exact(0,0)*axes(1,0)+du_exact(1,0)*axes(1,1);
		REAL parc1 = fabs(dx-dudx(0,0));
		REAL parc2 = fabs(dy-dudx(1,0));
		//Norma L2
		values[1] = pow(fabs(u[0] - u_exact[0]),(REAL)2.0);
		//seminorma
		values[2] = pow(parc1,(REAL)2.)+pow(parc2,(REAL)2.);
		//Norma Energia
		values[0] = values[1]+values[2];
		return;
	}
	//values[1] : eror em norma L2
	values[1]  = pow(sol[0] - u_exact[0],(REAL)2.0);
	//values[2] : erro em semi norma H1
	values[2]  = pow(dsol[0] - du_exact(0,0),(REAL)2.0);
	if(dudx.Rows()>1) values[2] += pow(dsol[1] - du_exact(1,0),(REAL)2.0);
	if(dudx.Rows()>2) values[2] += pow(dsol[2] - du_exact(2,0),(REAL)2.0);
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

void TPZMatOrthotropic::Normalize(TPZFMatrix &naxes){
	
	/** 
	 * os eixos devem vir por linhas e 
	 * orientados pela regra da m� direita
	 */
	
	int i;
	/**normalizando os eixos*/
	REAL norm2 = naxes(0,0)*naxes(0,0)+naxes(0,1)*naxes(0,1)+naxes(0,2)*naxes(0,2);
	if(norm2 == 0.){PZError << "TPZMatOrthotropic::Normalize: Eixo nulo nao e valido (eixo 1)\n"; exit(-1);}
	for(i=0;i<3;i++) naxes(0,i) /= sqrt(norm2);
	norm2 = naxes(1,0)*naxes(1,0)+naxes(1,1)*naxes(1,1)+naxes(1,2)*naxes(1,2);
	if(norm2 == 0.){PZError << "TPZMatOrthotropic::Normalize: Eixo nulo nao e valido (eixo 2)\n"; exit(-1);}
	for(i=0;i<3;i++) naxes(1,i) /= sqrt(norm2);
	norm2 = naxes(2,0)*naxes(2,0)+naxes(2,1)*naxes(2,1)+naxes(2,2)*naxes(2,2);
	if(norm2 != 0.) for(i=0;i<3;i++) naxes(2,i) /= sqrt(norm2);
	/**verificando a ortogonalidade dos dois primeiros eixos pelo produto vetorial*/
	REAL componente_K =  naxes(0,0)*naxes(1,1) - naxes(0,1)*naxes(1,0);
	REAL componente_J = -naxes(0,0)*naxes(1,2) + naxes(0,2)*naxes(1,0);
	REAL componente_I =  naxes(0,1)*naxes(1,2) - naxes(0,2)*naxes(1,1);
	/**primeiro teste*/
	if(componente_I == 0. && componente_J == 0. && componente_K == 0.){
		PZError << "TPZMatOrthotropic::Normalize: Os dois primeiros eixos nao devem ser paralelos\n";
		PZError << "Programa abortado\n";
		exit(-1);
	}
	/**aqui os dois primeiros vetores s� ortogonais*/
	/**o terceiro eixo �ortogonal aos dois primeiros*/
	norm2 = componente_I*componente_I+componente_J*componente_J+componente_K*componente_K;
	norm2 = sqrt(norm2);
	naxes(2,0) = componente_I/norm2;
	naxes(2,1) = componente_J/norm2;
	naxes(2,2) = componente_K/norm2;
	/**caso o segundo eixo n� seja ortogonal ao primeiro nem ao terceiro*/
	naxes(1,2) =  naxes(2,0)*naxes(0,1) - naxes(2,1)*naxes(0,0);//K
	naxes(1,1) = -naxes(2,0)*naxes(0,2) + naxes(2,2)*naxes(0,0);//J
	naxes(1,0) =  naxes(2,1)*naxes(0,2) - naxes(2,2)*naxes(0,1);//I
	naxes.Print(" * * * Eixos das fibras * * *");
}

//REAL norm = pow(norm2,0.5);
//REAL norm = sqrt(double(norm2));
//cout << "\n\nnorm = " << norm << endl;






/*
 REAL dphix = axes(0,0)*dsol[0]+axes(1,0)*dsol[1]+axes(2,0)*dsol[2];
 REAL dphiy = axes(0,1)*dsol[0]+axes(1,1)*dsol[1]+axes(2,1)*dsol[2];
 REAL dphiz = axes(0,2)*dsol[0]+axes(1,2)*dsol[1]+axes(2,2)*dsol[2];
 
 //values[2]  = pow(dphix - du_exact(0,0),2.0);
 //values[2] += pow(dphiy - du_exact(1,0),2.0);
 //values[2] += pow(dphiz - du_exact(2,0),2.0);
 if(dudx.Rows()<3) {
 REAL dx = du_exact(0,0)*axes(0,0);//+du_exact(1,0)*axes(0,1);
 REAL dy;
 if(dudx.Rows()>1)
 dx += du_exact(1,0)*axes(0,1);
 dy = du_exact(0,0)*axes(1,0)+du_exact(1,0)*axes(1,1);
 //values[1] : eror em norma L2
 values[1]  = pow(sol[0] - u_exact[0],2.0);
 //values[2] : erro em semi norma H1
 values[2]  = pow(dx - dudx(0,0),2.0);
 if(dudx.Rows()>1) values[2] += pow(dy - dudx(1,0),2.0);
 //values[0] : erro em norma H1 <=> norma Energia
 values[0]  = values[1]+values[2];
 return;
 }
 */

/*
 // dudx - dudy
 DSolxy[0][0] = DSol(0,0)*axes(0,0)+DSol(1,0)*axes(0,1)+DSol(2,0)*axes(0,2);
 DSolxy[1][0] = DSol(0,0)*axes(1,0)+DSol(1,0)*axes(1,1)+DSol(2,0)*axes(1,2);
 DSolxy[2][0] = DSol(0,0)*axes(2,0)+DSol(1,0)*axes(2,1)+DSol(2,0)*axes(2,2);
 // dvdx - dvdy
 DSolxy[0][1] = DSol(0,1)*axes(0,0)+DSol(1,1)*axes(0,1)+DSol(2,1)*axes(0,2);
 DSolxy[1][1] = DSol(0,1)*axes(1,0)+DSol(1,1)*axes(1,1)+DSol(2,1)*axes(1,2);
 DSolxy[2][1] = DSol(0,1)*axes(2,0)+DSol(1,1)*axes(2,1)+DSol(2,1)*axes(2,2);
 //
 DSolxy[0][2] = DSol(0,2)*axes(0,0)+DSol(1,2)*axes(0,1)+DSol(2,2)*axes(0,2);
 DSolxy[1][2] = DSol(0,2)*axes(1,0)+DSol(1,2)*axes(1,1)+DSol(2,2)*axes(1,2);
 DSolxy[2][2] = DSol(0,2)*axes(2,0)+DSol(1,2)*axes(2,1)+DSol(2,2)*axes(2,2);
 */
