/**
 * @file
 * @brief Contains implementations of the TPZMatHyperElastic methods.
 */

#include "pzmathyperelastic.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>

#include <cmath>
using namespace std;

TPZMatHyperElastic::TPZMatHyperElastic(int nummat,STATE e,STATE mu,STATE nu,
									   STATE lambda,STATE coef1,STATE coef2,STATE coef3) : 
TPZRegisterClassId(&TPZMatHyperElastic::ClassId),
TPZMaterial(nummat)
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

int TPZMatHyperElastic::NStateVariables() const {
	return 3;//3 deslocamentos
}

void TPZMatHyperElastic::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	TPZMaterial::Print(out);
}

void TPZMatHyperElastic::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZManVector<REAL,3> &x = data.x;
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	TPZFMatrix<STATE> &dsol=data.dsol[0];
	
	if(fForcingFunction) {
		TPZManVector<STATE> res(3);
		fForcingFunction->Execute(x,res);
		fXf[0] = res[0];
		fXf[1] = res[1];
		fXf[2] = res[2];
	}
	TFad<9, TFad<9,STATE> > U;
	ComputeEnergy(fLambda,fNu,dsol,U);
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
/*    for(int ii=0;ii<3;ii++)
        for(int jj=0;jj<3;jj++) {
            fK2[ii][jj] = (STATE)0.;
            fK3[ii][jj] = (STATE)0.;
            fK4[ii][jj] = (STATE)0.;
            fK6[ii][jj] = (STATE)0.;
            fK7[ii][jj] = (STATE)0.;
            fK8[ii][jj] = (STATE)0.;
        }
*/
}

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

void TPZMatHyperElastic::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
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
			
			TPZFMatrix<STATE> F(DSol),Ft(3,3,0.0),I(3,3,0.),S(3,3,0.);
			F(0,0)+=1;
			F(1,1)+=1;
			F(2,2)+=1;//a diagonal e' sempre a mesma
			F.Transpose(&Ft);
			TPZFMatrix<STATE> B = F*Ft;
			
			STATE ux,uy,uz,vx,vy,vz,wx,wy,wz;
			ux = DSol(0,0);
			uy = DSol(1,0);
			uz = DSol(2,0);
			
			vx = DSol(0,1);
			vy = DSol(1,1);
			vz = DSol(2,1);
			
			wx = DSol(0,2);
			wy = DSol(1,2);
			wz = DSol(2,2);
			
			STATE J = (ux+1.)*(vy+1.)*(wz+1.) + vx*wy*uz + wx*uy*vz - wx*(vy+1.)*uz - (ux+1.)*wy*vz - vx*uy*(wz+1.);
			
			I(0,0)=1.0;
			I(1,1)=1.0;
			I(2,2)=1.0;
			TPZFMatrix<STATE> sigmaF = (STATE(fNu/J)*B+STATE((fLambda*STATE(0.5)*(J*J-1.0)-fNu)/J)*I);
			STATE trsigma = sigmaF(0,0)+ sigmaF(1,1)+ sigmaF(2,2);
			S = sigmaF - STATE(trsigma/3.0)*I;
			int i,j;
			STATE J2=0.0;
			for(i=0;i<3;i++)  for(j=0;j<3;j++) J2 += S(i,j)* S(i,j);
			Solout[0] = sqrt(3.0*J2);
			
			
		}
		else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}

void TPZMatHyperElastic::Errors(TPZVec<REAL> &/*x*/,TPZVec<STATE> &u,
								TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
								TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
	//TPZVec<REAL> sol(1),dsol(3);
	TPZManVector<STATE> sol(3),dsol(9);
	Solution(u,dudx,axes,2,sol);
	Solution(u,dudx,axes,3,dsol);
	//values[1] : erro em norma L2
	values[1]  = pow(sol[0] - u_exact[0],(STATE)2.0);
	values[1] += pow(sol[1] - u_exact[1],(STATE)2.0);
	values[1] += pow(sol[2] - u_exact[2],(STATE)2.0);
	//values[2] : erro em semi norma H1
	int k=0;
	values[2] = 0.;
	for(int i=0;i<3;i++) {
		values[2] += pow(dsol[k++] - du_exact(0,i),(STATE)2.0);
		values[2] += pow(dsol[k++] - du_exact(1,i),(STATE)2.0);
		values[2] += pow(dsol[k++] - du_exact(2,i),(STATE)2.0);
	}
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

void TPZMatHyperElastic::ContributeBC(TPZMaterialData &data,
                                   REAL weight,
                                   TPZFMatrix<STATE> &ek,
                                   TPZFMatrix<STATE> &ef,
                                   TPZBndCond &bc){
	TPZFMatrix<REAL> &phi = data.phi;
	
	const STATE BIGNUMBER  = 1.e12;
	
	const int phr = phi.Rows();
	int in,jn,idf,jdf;
	STATE v2[3];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	v2[2] = bc.Val2()(2,0);
	TPZFMatrix<STATE> &v1 = bc.Val1();
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
	switch (bc.Type()) {
		case 0: // Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(3*in+0,0) += BIGNUMBER * v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += BIGNUMBER * v2[1] * phi(in,0) * weight;        
				ef(3*in+2,0) += BIGNUMBER * v2[2] * phi(in,0) * weight;        
				
				for (jn = 0 ; jn < phr; jn++) {
					ek(3*in+0,3*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
					ek(3*in+1,3*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
					ek(3*in+2,3*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight;
				}//jn
			}//in
			break;
			
		case 1: // Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
			}//in
			break;
		case 2: // Mixed condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
				for(jn=0; jn<phi.Rows(); jn++)
				{
					for(idf=0; idf<3; idf++) for(jdf=0; jdf<3; jdf++)
					{
						ek(3*in+idf,3*jn+jdf) += bc.Val1()(idf,jdf);
					}
				}
			}//in
			break;
		case 3: // Directional Null Dirichlet - displacement is set to null in the non-null vector component direction
			for(in = 0 ; in < phr; in++) {
				ef(3*in+0,0) += BIGNUMBER * (0. - data.sol[0][0]) * v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += BIGNUMBER * (0. - data.sol[0][1]) * v2[1] * phi(in,0) * weight;        
				ef(3*in+2,0) += BIGNUMBER * (0. - data.sol[0][2]) * v2[2] * phi(in,0) * weight;        
				for (jn = 0 ; jn < phr; jn++) {
					ek(3*in+0,3*jn+0) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[0];
					ek(3*in+1,3*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[1];
					ek(3*in+2,3*jn+2) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2[2];
				}//jn
			}//in
			break;
			
		case 4: // stressField Neumann condition
			for(in = 0; in < 3; in ++)
				v2[in] = - ( v1(in,0) * data.normal[0] +
							v1(in,1) * data.normal[1] +
							v1(in,2) * data.normal[2] );
			// The normal vector points towards the neighbour. The negative sign is there to 
			// reflect the outward normal vector.
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(3*in+0,0) += v2[0] * phi(in,0) * weight;
				ef(3*in+1,0) += v2[1] * phi(in,0) * weight;
				ef(3*in+2,0) += v2[2] * phi(in,0) * weight;
				//cout << "normal:" << data.normal[0] << ' ' << data.normal[1] << ' ' << data.normal[2] << endl;
				//cout << "val2:  " << v2[0]          << ' ' << v2[1]          << ' ' << v2[2]          << endl;
			}
			break;
		default:
			PZError << "TPZElastitity3D::ContributeBC error - Wrong boundary condition type" << std::endl;
	}//switch
	
}//method

//
//void TPZMatHyperElastic::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc){
//	
//	TPZFMatrix<REAL> &phi = data.phi;
//    int numbersol = data.sol.size();
//    if (numbersol != 1) {
//        DebugStop();
//    }
//	TPZVec<STATE> &sol=data.sol[0];
//	
//	if(bc.Material() != this){
//		PZError << "TPZMatHyperElastic.ContributeBC : this material don't exists \n";
//	}
//	
//	if(bc.Type() < 0 && bc.Type() > 2){
//		PZError << "ContributeBC.aplybc, unknown boundary condition type : "<<bc.Type() << endl;
//	}
//	
//	int ndof = NStateVariables();
//	int nnod = ek.Rows()/ndof;
//	int r = ndof;
//	
//	int idf,jdf,in,jn;
//	switch(bc.Type()){
//		case 0:
//			for(in=0 ; in<nnod ; ++in){
//				for(idf = 0;idf<r;idf++) {
//					(ef)(in*r+idf,0) += gBigNumber*phi(in,0)*(bc.Val2()(idf,0)-sol[idf])*weight;
//				}
//				for(jn=0 ; jn<nnod ; ++jn) {
//					for(idf = 0;idf<r;idf++) {
//						ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
//					}
//				}
//			}
//			break;
//			
//		case 1:
//			for(in=0 ; in<nnod ; ++in){
//				for(idf = 0;idf<r;idf++) {
//					//(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0)-sol[idf]);
//					(ef)(in*r+idf,0) += weight*phi(in,0)*(bc.Val2()(idf,0));
//				}
//			}
//			break;
//			
//		case 2:
//			for(in=0 ; in<nnod ; ++in){
//				for(idf = 0;idf<r;idf++) {
//					for (jdf=0; jdf<r; jdf++){
//						(ef)(in*r+idf,0) += phi(in,0)*bc.Val1()(idf,jdf)*(bc.Val2()(jdf,0)-sol[jdf])*weight;
//					}
//					for(jn=0 ; jn<nnod ; ++jn) {
//						for(idf = 0;idf<r;idf++) {
//							for(jdf = 0;jdf<r;jdf++) {
//								ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
//							}
//						}
//					}
//				}
//				
//            }
//			
//	}//fim switch
//}


/** The function below makes the correspondence between the dsol vector and a matrix ordered F operator */
inline int ith(const int i, const int j)
{
	return i*3+j;
}

void TPZMatHyperElastic::ContributeEnergy(TPZVec<REAL> &x,
										  TPZVec<FADFADREAL> &sol, TPZVec<FADFADREAL> &dsol,
										  FADFADREAL &U, REAL weight)
{
    DebugStop();
    /*
	FADFADREAL J, TrC; // J = det(F); TrC = Trace(C)
	FADFADREAL DiagF0(dsol[    0]);
	DiagF0.val().val() += 1.;// element [0][0]
	
	FADFADREAL DiagF1(dsol[  3+1]);
	DiagF1.val().val() += 1.;// element [0][0]
	
	FADFADREAL DiagF2(dsol[2*3+2]);
	DiagF2.val().val() += 1.;// element [0][0]
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
 */
}

void TPZMatHyperElastic::ContributeBCEnergy(TPZVec<REAL> & x,
											TPZVec<FADFADREAL> & sol, FADFADREAL &U,
											REAL weight, TPZBndCond &bc)
{
    DebugStop();
    /*
	if(bc.Material() != this){
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
     */
}

void TPZMatHyperElastic::ComputeEnergy(STATE lambda, STATE mu,  TPZFMatrix<STATE> &dsol, TFad<9,TFad<9,STATE> > &energy) {
    DebugStop();
    /*
	TFad<9,TFad<9,STATE> > tensor[3][3],J,TrC;
	int i,j;
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
	
	J = tensor[0][0] * tensor[1][1] * tensor[2][2] + 
    tensor[0][1] * tensor[1][2] * tensor[2][0] +
    tensor[0][2] * tensor[1][0] * tensor[2][1] -
    tensor[0][2] * tensor[1][1] * tensor[2][0] -
    tensor[0][1] * tensor[1][0] * tensor[2][2] -
    tensor[0][0] * tensor[1][2] * tensor[2][1]; //  J = det(F)
	
	energy = (J*J - TFad<9,STATE>(1.)) * TFad<9,STATE>(lambda/4.) -
    log( J ) * TFad<9,STATE>((lambda/2.+mu)) +
    (TrC - TFad<9,STATE>(3.)) * TFad<9,STATE>(mu/2.);
	*/
}


int TPZMatHyperElastic::ClassId() const{
    return Hash("TPZMatHyperElastic") ^ TPZMaterial::ClassId() << 1;
}
