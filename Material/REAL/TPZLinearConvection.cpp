/**
 * \file
 * @brief Contains implementations of the TPZLinearConvection methods.
 */

#include "TPZLinearConvection.h"
#include "pzbndcond.h"
#include "pzerror.h"

using namespace std;

TPZLinearConvection::TPZLinearConvection(TPZLinearConvection & copy) :
TPZRegisterClassId(&TPZLinearConvection::ClassId), TPZMaterial(copy){
    fConvect[0] = copy.fConvect[0];
    fConvect[1] = copy.fConvect[1];
}
TPZLinearConvection::TPZLinearConvection(int id, TPZVec<STATE>& conv) :
TPZRegisterClassId(&TPZLinearConvection::ClassId), TPZMaterial(id){
    fConvect[0] = conv[0];
    fConvect[1] = conv[1];
}

void TPZLinearConvection::SetData(istream &data) {
	PZError << "TPZMaterial::SetData is called.\n";
	data >> fConvect[0] >> fConvect[1];
}
TPZMaterial * TPZLinearConvection::NewMaterial() {
	TPZLinearConvection *result = new TPZLinearConvection(*this);
	
	return result;
}
void TPZLinearConvection::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,
								   TPZVec<STATE> &Solout){
	if(var == 1) {
		Solout.Resize(2);
		Solout[0] = Sol[0]*fConvect[0];
		Solout[1] = Sol[0]*fConvect[1];
	} else TPZMaterial::Solution(Sol,DSol,axes,var,Solout);
}
int TPZLinearConvection::NSolutionVariables(int index) {
	if(index == 1) return 2;
	return TPZMaterial::NSolutionVariables(index);
}
int TPZLinearConvection::VariableIndex(const std::string &name) {
	if(!strcmp(name.c_str(),"flux")) return 1;
	return TPZMaterial::VariableIndex(name);
}
void TPZLinearConvection::Print(std::ostream & out) {
    TPZMaterial::Print(out);
    out << "Convection : " << fConvect[0] << ' ' << fConvect[1] << endl;
}
int TPZLinearConvection::NStateVariables()  const {
    return 1;
}
int TPZLinearConvection::Dimension() const {
    return 2;
}
void TPZLinearConvection::Contribute(TPZMaterialData &data, REAL weight,
									 TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZFMatrix<REAL> &phi = data.phi;
	TPZFMatrix<REAL> &daxesdksi=data.jacinv;
	TPZFMatrix<REAL> &axes=data.axes;
	
    STATE convectax[2];
    convectax[0] = fConvect[0]*axes(0,0)+fConvect[1]*axes(0,1);
    convectax[1] = fConvect[0]*axes(1,0)+fConvect[1]*axes(1,1);
    STATE sconvect[2] = {1.,1.};
    if(convectax[0] < 0.) sconvect[0] = -1.;
    if(convectax[1] < 0.) sconvect[1] = -1.;
	//    REAL convectnorm = sqrt(convectax[0]*convectax[0]+convectax[1]*convectax[1]);
    STATE delax[2];
    delax[0] = sqrt(daxesdksi(0,0)*daxesdksi(0,0)+daxesdksi(1,0)*daxesdksi(1,0));
    delax[1] = sqrt(daxesdksi(0,1)*daxesdksi(0,1)+daxesdksi(1,1)*daxesdksi(1,1));
    int nshape = phi.Rows();
    int in,jn;
    for(in=0; in<nshape; in++) {
        for(jn=0; jn<nshape; jn++) {
            ek(in,jn) += weight*(
								 delax[0]*sconvect[0]*dphi(0,in)*(convectax[0]*dphi(0,jn)+convectax[1]*dphi(1,jn))
								 +delax[1]*sconvect[1]*dphi(1,in)*(convectax[0]*dphi(0,jn)+convectax[1]*dphi(1,jn))
								 //delax[0]/convectnorm*(dphi(0,in)*convectax[0]+dphi(1,in)*convectax[1])*(dphi(0,jn)*convectax[0]+dphi(1,jn)*convectax[1])
								 -dphi(0,in)*convectax[0]*phi(jn,0)
								 -dphi(1,in)*convectax[1]*phi(jn,0)
								 );
        }
    }
}

void TPZLinearConvection::ContributeBC(TPZMaterialData &data,REAL weight,
									   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	TPZFMatrix<REAL> &phi = data.phi;
	
	if(bc.Material() != this){
		PZError << "TPZMat1dLin.apply_bc warning : this material didn't create the boundary condition!\n";
	}
	
	if(bc.Type() < 0 && bc.Type() > 2){
		PZError << "TPZMat1dLin.aplybc, unknown boundary condition type :"  <<
		bc.Type() << " boundary condition ignored\n";
	}
	
	int numdof = NStateVariables();
	int numnod = ek.Rows()/numdof;
	int r = numdof;
	
	int idf,jdf,in,jn;
	switch(bc.Type()){
			
		case 0:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += gBigNumber*phi(in,0)*bc.Val2()(idf,0)*weight;
				}
				for(jn=0 ; jn<numnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						ek(in*r+idf,jn*r+idf) += gBigNumber*phi(in,0)*phi(jn,0)*weight;
					}
				}
			}
			break;
			
		case 1:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
				}
			}
			break;
			
		case 2:
			for(in=0 ; in<numnod ; ++in){
				for(idf = 0;idf<r;idf++) {
					(ef)(in*r+idf,0) += phi(in,0)*bc.Val2()(idf,0)*weight;
				}
				for(jn=0 ; jn<numnod ; ++jn) {
					for(idf = 0;idf<r;idf++) {
						for(jdf = 0;jdf<r;jdf++) {
							ek(in*r+idf,jn*r+jdf) += bc.Val1()(idf,jdf)*phi(in,0)*phi(jn,0)*weight;
						}
					}
				}
			}//fim switch
	}
}

int TPZLinearConvection::ClassId() const{
    return Hash("TPZLinearConvection") ^ TPZMaterial::ClassId() << 1;
}
