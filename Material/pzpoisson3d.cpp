/**
 * @file
 * @brief Contains implementations of the TPZMatPoisson3d methods.
 */

#include "pzpoisson3d.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"

#include <cmath>

#ifdef PZ_LOG
static TPZLogger logger("pz.material.poisson3d");
#endif


using namespace std;
STATE TPZMatPoisson3d::gAlfa = 0.5;

TPZMatPoisson3d::TPZMatPoisson3d(int nummat, int dim) : TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
TPZMaterial(nummat), fXf(0.), fDim(dim), fSD(0.) {
    if(dim < 1)
    {
        DebugStop();
    }
    fTensorPerm ={{1,0,0},{0,1,0},{0,0,1}};
    fInvPerm ={{1,0,0},{0,1,0},{0,0,1}};
	fC = 0.;
	fConvDir[0] = 0.;
	fConvDir[1] = 0.;
	fConvDir[2] = 0.;
	fPenaltyConstant = 1000.;
	this->SetNonSymmetric();
	this->SetNoPenalty();
  fNeumann=false;
}

TPZMatPoisson3d::TPZMatPoisson3d():TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
TPZMaterial(), fXf(0.), fDim(1), fSD(0.){
	fC = 0.;
    fTensorPerm ={{1,0,0},{0,1,0},{0,0,1}};
    fInvPerm ={{1,0,0},{0,1,0},{0,0,1}};
	fConvDir[0] = 0.;
	fConvDir[1] = 0.;
	fConvDir[2] = 0.;
	fPenaltyConstant = 1000.;
	this->SetNonSymmetric();
	this->SetNoPenalty();
  fNeumann=false;
}

TPZMatPoisson3d::TPZMatPoisson3d(const TPZMatPoisson3d &copy):TPZRegisterClassId(&TPZMatPoisson3d::ClassId),
TPZMaterial(copy){
	this->operator =(copy);
}

TPZMatPoisson3d & TPZMatPoisson3d::operator=(const TPZMatPoisson3d &copy){
	TPZMaterial::operator = (copy);
	fXf  = copy.fXf;
	fDim = copy.fDim;
	fTensorPerm = copy.fTensorPerm;
	fInvPerm = copy.fInvPerm;
	fC   = copy.fC;
	for (int i = 0; i < 3; i++) fConvDir[i] = copy.fConvDir[i];
	fSymmetry = copy.fSymmetry;
	fSD = copy.fSD;
	fPenaltyConstant = copy.fPenaltyConstant;
	this->fPenaltyType = copy.fPenaltyType;
	return *this;
}

void TPZMatPoisson3d::SetParameters(STATE diff, REAL conv, TPZVec<REAL> &convdir) {
	fTensorPerm.Resize(3,3);
	fTensorPerm.Zero();
    fTensorPerm(0,0) = fTensorPerm(1,1) = fTensorPerm(2,2) = diff;
	fC = conv;
	int d;
	for(d=0; d<fDim; d++) fConvDir[d] = convdir[d];
}

void TPZMatPoisson3d::SetParameters(TPZFNMatrix<9,STATE> tensorPerm,TPZFNMatrix<9,STATE> invPerm,REAL conv, TPZVec<REAL> &convdir){
    if(tensorPerm.Rows() != 3 || tensorPerm.Cols() != 3 || invPerm.Rows() != 3 || invPerm.Cols() != 3) {
        std::cout << "ERROR\nvoid SetParameters(TPZFMatrix<STATE> tensorDiff,TPZFMatrix<STATE> invPerm,REAL conv, TPZVec<REAL> &convdir): Invalid matrix dimension\n";
    }
    fTensorPerm = tensorPerm;
    fInvPerm = invPerm;
    fC = conv;
    int d;
    for(d=0; d<fDim; d++) fConvDir[d] = convdir[d];
}

void TPZMatPoisson3d::GetParameters(TPZFNMatrix<9,STATE> &tensorPerm, TPZFNMatrix<9,STATE> &invPerm, REAL &conv, TPZVec<REAL> &convdir) {
	tensorPerm = fTensorPerm;
	invPerm = fInvPerm;
	conv = fC;
	int d;
	for(d=0; d<fDim; d++) convdir[d] = fConvDir[d];
}

void TPZMatPoisson3d::SetPermeability(STATE perm){
    fTensorPerm.Resize(3,3); fTensorPerm.Zero();
    fInvPerm.Resize(3,3); fInvPerm.Zero();
    for (int i=0; i<3; i++) {
        fTensorPerm(i,i) = perm;
        fInvPerm(i,i) = 1./perm;
    }
}

void TPZMatPoisson3d::SetPermeabilityTensor(const TPZFNMatrix<9,STATE> K, const TPZFNMatrix<9,STATE> invK){

    if(K.Rows() != 3 || invK.Rows() != 3)
    {
        if(K.Rows() != 2 || invK.Rows() != 2) {
            std::cout << "ERROR: Insert a 3x3 or a 2x2 permeability tensor";
            DebugStop();
        }
        for(int iind =0; iind < 2 ; iind++) for(int jind = 0; jind <2 ; jind++){
            fTensorPerm(iind,jind) = K.GetVal(iind,jind);
            fInvPerm(iind,jind) = invK.GetVal(iind,jind);
        }
    }
    else {
        fTensorPerm = K;
        fInvPerm = invK;
    }
}

void TPZMatPoisson3d::GetPermeability(TPZFNMatrix<9,STATE> &K){
    K = fTensorPerm;
}

void TPZMatPoisson3d::GetInvPermeability(TPZFNMatrix<9,STATE> &invK){
    invK = fInvPerm;
}

/// return the permeability and compute it if there is permeability function
void TPZMatPoisson3d::GetPermeabilities(TPZVec<REAL> &x, TPZFNMatrix<9,STATE> &PermTensor, TPZFNMatrix<9,STATE> &InvPermTensor)
{
    this->GetPermeability(PermTensor);
    this->GetInvPermeability(InvPermTensor);
}

void TPZMatPoisson3d::GetParameters(STATE &diff, REAL &conv, TPZVec<REAL> &convdir) {
#ifdef PZDEBUG
    std::cout <<"WARNING: TPZMatPoisson3d::GetParameters: Permeability is a tensor. This method only works if one is working with an identity permeability tensor. If that is not one's case, the result is wrong.";
#endif
    diff = AVGK();
    conv = fC;
    int d;
    for(d=0; d<fDim; d++) convdir[d] = fConvDir[d];
}

TPZMatPoisson3d::~TPZMatPoisson3d() {
}

//int TPZMatPoisson3d::NStateVariables() const {
//    return 1;
//}

void TPZMatPoisson3d::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
    out << "Permeability Tensor "; fTensorPerm.Print(out);
	out << "Convection coeficient fC " << fC << endl;
	out << "Convection direction " << fConvDir[0] << ' ' << fConvDir[1] << ' ' <<  fConvDir[2] << endl;
	out << "Forcing vector fXf " << fXf << endl;
	out << "Penalty constant " << fPenaltyConstant << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMatPoisson3d::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
    
    if(fNeumann){
        
        LocalNeumanContribute(data , weight , ek, ef);
        
        return;
    }

    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();

    STATE fXfLoc = fXf;

    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fXfLoc = -res[0];
    }
    REAL delx = 0.;
    STATE ConvDirAx[3] = {0.};
    if(fC != 0.0) {
        int di,dj;
        delx = 0.;
        for(di=0; di<fDim; di++) {
            for(dj=0; dj<fDim; dj++) {
                delx = (delx<fabs(jacinv(di,dj))) ? fabs(jacinv(di,dj)) : delx;
            }
        }
        delx = 2./delx;
        
        
        switch(fDim) {
            case 1:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                break;
            case 2:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
                break;
            case 3:
                ConvDirAx[0] = axes(0,0)*fConvDir[0]+axes(0,1)*fConvDir[1]+axes(0,2)*fConvDir[2];
                ConvDirAx[1] = axes(1,0)*fConvDir[0]+axes(1,1)*fConvDir[1]+axes(1,2)*fConvDir[2];
                ConvDirAx[2] = axes(2,0)*fConvDir[0]+axes(2,1)*fConvDir[1]+axes(2,2)*fConvDir[2];
                break;
            default:
                PZError << "TPZMatPoisson3d::Contribute dimension error " << fDim << endl;
        }
    }

    //Kdphi = K*dphi
    TPZFMatrix<STATE> Kdphi(fDim,phr,0);
    for(int in = 0 ; in < fDim ; in++) for(int jn = 0; jn < phr; jn++) for(int kn =0; kn <fDim; kn++)
        Kdphi(in,jn) += fTensorPerm.GetVal(in,kn)*dphi(kn,jn);

    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        STATE dphiic = 0;
        for(kd = 0; kd<fDim; kd++) dphiic += ConvDirAx[kd]*(STATE)dphi(kd,in);
        ef(in, 0) += - (STATE)weight * fXfLoc * ( (STATE)phi(in,0) + (STATE)(0.5*delx*fC)*fSD*dphiic );
        for( int jn = 0; jn < phr; jn++ ) {

            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight * (
                                       + (STATE)( dphi(kd,in) * Kdphi(kd,jn) )
                                       - (STATE)(fC* dphi(kd,in) * phi(jn)) * ConvDirAx[kd]
                                       + (STATE)(0.5 * delx * fC * dphi(kd,jn)) * fSD * dphiic * ConvDirAx[kd]
                                       );
            }
        }
    }

    if (this->IsSymetric()){    
        if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }


}

//
void TPZMatPoisson3d::LocalNeumanContribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    /** monta a matriz
     |Sk Ck^T |  = |bk|
     |Ck  0   |    |uk|
     
     **/
    int newRows=ek.Rows()+1;
    int newCols=ek.Cols()+1;
    ek.Resize( newRows, newCols );
    ef.Resize(newRows,1);
   
    STATE fXfLoc = fXf;
    if(fForcingFunction) {                           // phi(in, 0) = phi_in
        TPZManVector<STATE> res(1);
        fForcingFunction->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
        fXfLoc = res[0];
    }
    
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    TPZFMatrix<REAL> &axes = data.axes;
    TPZFMatrix<REAL> &jacinv = data.jacinv;
    TPZVec<STATE> &sol=data.sol[0];
    TPZGradSolVec &dsol=data.dsol;

    
    int phr = phi.Rows();
    int nlinhaek=ek.Rows();
    int ncolek=ek.Cols();

    //Kdphi = K*dphi
    TPZFMatrix<STATE> Kdphi(3,3,0);
    for(int in = 0 ; in < 3 ; in++) for(int jn = 0; jn < 3; jn++) for(int kn =0; kn <3; kn++)
                Kdphi(in,jn) += fTensorPerm.GetVal(in,kn)*dphi(kn,jn);
    
    int i,j;
    for( int in = 0; in < phr; in++ ){
        for( int jn = 0; jn < phr; jn++ ) {
            for(int kd=0; kd<fDim; kd++) {
            ek(in,jn) += (STATE)weight * ((STATE)( dphi(kd,in) * Kdphi(kd,jn) ));//gradphi_i.gradphi_j
            }
        }
        ek(in,ncolek-1) += (STATE)weight * phi(in);
        ek(nlinhaek-1,in) += (STATE)weight * phi(in);
    }
    
    for( int in = 0; in < phr; in++ ) {
        for( int jn = 0; jn < phr; jn++ ) {
            for(int kd=0; kd<fDim; kd++) {
                ef(in,0) += (STATE)weight * ( dphi(kd,in) * dsol[0](kd,jn) );
            }
        }
    }
 ef(phr,0)+=(STATE)weight * (sol[0]);
    
}


///

void TPZMatPoisson3d::ContributeBC(TPZMaterialData &data,REAL weight,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &axes = data.axes;
	int phr = phi.Rows();
	short in,jn;
	STATE v2[1];
	v2[0] = bc.Val2()(0,0);
	
	if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
		TPZManVector<STATE,1> res(1);
		bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
//        if(fabs(res[0]) > 1.e-6)
//        {
//            bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
//            DebugStop();
//        }
		v2[0] = res[0];
	}

	switch (bc.Type()) {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phr; in++) {
				ef(in,0) += (STATE)(gBigNumber* phi(in,0) * weight) * v2[0];
				for (jn = 0 ; jn < phr; jn++) {
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		case 1 :			// Neumann condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in,0) += v2[0] * (STATE)(phi(in,0) * weight);
			}
			break;
		case 2 :		// mixed condition
			for(in = 0 ; in < phi.Rows(); in++) {
				ef(in, 0) += v2[0] * (STATE)(phi(in, 0) * weight);
				for (jn = 0 ; jn < phi.Rows(); jn++) {
					ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi(in,0) * phi(jn,0) * weight);     // peso de contorno => integral de contorno
				}
			}
			break;
		case 3: // outflow condition
			int id, il, jl;
			REAL normal[3];
			if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
			if (fDim == 2){
				normal[0] = axes(0,1);
				//normal[1] = axes(1,1);
			}
			if (fDim == 3){
				normal[0] = axes(0,2);
				normal[1] = axes(1,2);
				normal[2] = axes(2,2);
			}
			REAL ConvNormal = 0.;    
			for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];  
			if(ConvNormal > 0.) {
				for(il=0; il<phr; il++) {
					for(jl=0; jl<phr; jl++) {
						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
					}
				}
			}
			else{
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
			break;
	}
	
	if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
		if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}

/** Returns the variable index associated with the name */
int TPZMatPoisson3d::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
	if(!strcmp("Derivative",name.c_str()))      return  2;
	if(!strcmp("KDuDx",name.c_str()))           return  3;
	if(!strcmp("KDuDy",name.c_str()))           return  4;
	if(!strcmp("KDuDz",name.c_str()))           return  5;
	if(!strcmp("NormKDu",name.c_str()))         return  6;
	if(!strcmp("MinusKGradU",name.c_str()))     return  7;
	if(!strcmp("POrder",name.c_str()))          return  8;
	if(!strcmp("Laplac",name.c_str()))          return  9;
	if(!strcmp("Stress",name.c_str()))          return  10;    
	if(!strcmp("Flux",name.c_str()))            return  10;
	if(!strcmp("Pressure",name.c_str()))        return  11;
	
	if (!strcmp("ExactPressure", name.c_str()))   return  12;
	if(!strcmp("ExactSolution",name.c_str()))   return  12;
	if(!strcmp("ExactFlux",name.c_str()))       return  13;
	if(!strcmp("Divergence",name.c_str()))      return  14;
	if(!strcmp("ExactDiv",name.c_str()))        return  15;
	
	if(!strcmp("PressureOmega1",name.c_str()))  return  16;
	if(!strcmp("PressureOmega2",name.c_str()))  return  17;
	if(!strcmp("FluxOmega1",name.c_str()))      return  18;
    
    if(!strcmp("GradFluxX",name.c_str()))       return  19;
    if(!strcmp("GradFluxY",name.c_str()))       return  20;
     if(!strcmp("FluxL2",name.c_str()))            return  21;//Only To calculate l2 error
	 if(!strcmp("OrdemP",name.c_str()))        return  99;
	return TPZMaterial::VariableIndex(name);
}

int TPZMatPoisson3d::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return fDim;//arrumar o fluxo de hdiv para ser fdim tbem enquanto isso faco isso
	if ((var == 3) || (var == 4) || (var == 5) || (var == 6)) return 1;
	if (var == 7) return fDim;
	if (var == 8) return 1;
	if (var == 9) return 1;
	if (var==10) return fDim;
	if (var==11) return 1;
	
	if (var==12) return 1;
	if (var==13) return fDim;
	if (var==14) return 1;
	if (var==15) return 1;
	//teste de acoplamento
	if (var==16) return 1;
	if (var==17) return 1;
	if (var==18) return 3;
    if (var==19) return 3;
    if (var==20) return 3;
    if (var==21) return fDim;
	
	
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatPoisson3d::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
	
	TPZVec<STATE> pressure(1);
	TPZVec<REAL> pto(3);
	TPZFMatrix<STATE> flux(3,1);
	
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
	
   // Solout.Resize(this->NSolutionVariables(var));
    
#ifndef STATE_COMPLEX
	
	switch (var) {
	/*	case 7:
		{
			//			{ //MinusKGradU
			int id;
			TPZManVector<STATE> dsolp(3,0);
			dsolp[0] = data.dsol[0](0,0)*data.axes(0,0)+data.dsol[0](1,0)*data.axes(1,0);
			dsolp[1] = data.dsol[0](0,0)*data.axes(0,1)+data.dsol[0](1,0)*data.axes(1,1);			
			dsolp[2] = data.dsol[0](0,0)*data.axes(0,2)+data.dsol[0](1,0)*data.axes(1,2);
			for(id=0 ; id<fDim; id++) 
			{
				Solout[id] = -1. * this->fK * dsolp[id];
			}
		}
			break;*/
		case 8:
			Solout[0] = data.p;
			break;
		case 10:
			if (data.numberdualfunctions) {
				
				Solout[0]=data.sol[0][0];
				Solout[1]=data.sol[0][1];
                Solout[2]=data.sol[0][2];
				
			}
			else {
				this->Solution(data.sol[0], data.dsol[0], data.axes, 2, Solout);
			}
			
			break;
            
        case 21:
            for(int k=0;k<fDim;k++){
                Solout[k]=data.sol[0][k];
            }
			break;
            
		case 11:
			if (data.numberdualfunctions) {
				Solout[0]=data.sol[0][2];
			}
			else{
				Solout[0]=data.sol[0][0];
			}
			break;
			
		case 12:
				fExactSol->Execute(data.x,pressure,flux);
				
				Solout[0]=pressure[0];
			break;
		case 13:
				fExactSol->Execute(data.x,pressure,flux);
				
				Solout[0]=flux(0,0);
				Solout[1]=flux(1,0);
            break;
            
        case 14:
        {
			if (data.numberdualfunctions){
				Solout[0]=data.sol[0][data.sol[0].NElements()-1];
			}else{
                //Solout[0]=data.dsol[0](0,0)+data.dsol[0](1,1)+data.dsol[0](2,2);
                STATE val = 0.;
                for(int i=0; i<fDim; i++){
                    val += data.dsol[0](i,i);
                }
                Solout[0] = val;
            }
        }
            break;
          
        case 15:
        {
            fExactSol->Execute(data.x,pressure,flux);
            Solout[0]=flux(fDim,0);
        }
            break;

            
        case 16:
            if (data.numberdualfunctions) {
					Solout[0]=data.sol[0][2];
            }
            else {
                std::cout<<"Pressao somente em Omega1"<<std::endl;
                Solout[0]=0;//NULL;
            }
				
            break;
        
        case 17:
            if (!data.numberdualfunctions) {
                Solout[0]=data.sol[0][0];
            }
            else {
                std::cout<<"Pressao somente em omega2"<<std::endl;
                Solout[0]=0;//NULL;
            }
				
            break;
        case 18:
            if( data.numberdualfunctions){
                Solout[0]=data.sol[0][0];//fluxo de omega1
                Solout[1]=data.sol[0][1];
                //	Solout[2]=data.sol[2];
                return;
            }
        
        case 19:
            if(data.numberdualfunctions){
                Solout[0]=data.dsol[0](0,0);//fluxo de omega1
                Solout[1]=data.dsol[0](1,0);
                Solout[2]=data.dsol[0](2,0);
                return;
            }
        case 20:
            if( data.numberdualfunctions){
                Solout[0]=data.dsol[0](0,1);//fluxo de omega1
                Solout[1]=data.dsol[0](1,1);
                Solout[2]=data.dsol[0](2,1);
                return;
            }
            else {
                std::cout<<"Pressao somente em omega2"<<std::endl;
                Solout[0]=0;//NULL;
            }
            break;
        default:
           
            if (data.sol[0].size() == 4) {
                
                data.sol[0][0] = data.sol[0][2];
            }

            this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
            break;
    }
#endif
}

#include "pzaxestools.h"
void TPZMatPoisson3d::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
#ifndef STATE_COMPLEX
	Solout.Resize( this->NSolutionVariables( var ) );
	
	if(var == 1){
		Solout[0] = Sol[0];//function
		return;
	}
	if(var == 2) {
		int id;
		for(id=0 ; id<fDim; id++) {
			TPZFNMatrix<9,STATE> dsoldx;
			TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
			Solout[id] = dsoldx(id,0);//derivate
		}
		return;
	}//var == 2
	if (var == 3){ //KDuDx
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = 0;
		for(int index = 0; index < this->fTensorPerm.Cols(); index++)
		    Solout[0] += (this->fTensorPerm.GetVal(0,index))*dsoldx(index,0);
		return;
	}//var ==3
	if (var == 4){ //KDuDy
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
        Solout[0] = 0;
        for(int index = 0; index < this->fTensorPerm.Cols(); index++)
            Solout[0] += (this->fTensorPerm.GetVal(1,index))*dsoldx(index,0);
		return;
	}//var == 4 
	if (var == 5){ //KDuDz
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
        Solout[0] = 0;
        for(int index = 0; index < this->fTensorPerm.Cols(); index++)
            Solout[0] += (this->fTensorPerm.GetVal(2,index))*dsoldx(index,0);
		return;
	}//var == 5
	if (var == 6){ //NormKDu
		int id;
        TPZFNMatrix<9,STATE> dsoldx;
        TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		STATE val = 0., dsoldxx;
		for(id=0 ; id<fDim; id++){
		    dsoldxx = 0.;
            for(int index = 0; index < fDim; index++)
                dsoldxx += (this->fTensorPerm.GetVal(id,index))*dsoldx(index);
            val += dsoldxx*dsoldxx;
		}
		Solout[0] = sqrt(val);
		return;
	}//var == 6
	if (var == 7){ //MinusKGradU
		int id;
		//REAL val = 0.;
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
        for(id=0 ; id<fDim; id++) Solout[id] =0.;
		for(id=0 ; id<fDim; id++) {
            for(int index = 0; index < this->fTensorPerm.Cols(); index++)
                Solout[id] -= (this->fTensorPerm.GetVal(id,index))*dsoldx(index,0);
		}
		return;
	}//var == 7  
	if(var == 9){//Laplac
		Solout.Resize(1);
		Solout[0] = DSol(2,0);
		return;
	}//Laplac
	
#endif
	TPZMaterial::Solution(Sol, DSol, axes, var, Solout);
	
}//method

void TPZMatPoisson3d::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
							 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, 
							 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
	values.Resize(NEvalErrors());
    values.Fill(0.0);
    TPZVec<STATE> dudX(3,0);
	TPZManVector<STATE> dudxEF(1,0.), dudyEF(1,0.),dudzEF(1,0.);
	this->Solution(u,dudx,axes, this->VariableIndex("KDuDx"), dudxEF);
    this->Solution(u,dudx,axes, this->VariableIndex("KDuDy"), dudyEF);
    this->Solution(u,dudx,axes, this->VariableIndex("KDuDz"), dudzEF);
    dudX[0] = dudxEF[0]; dudX[1] = dudyEF[0]; dudX[2] = dudzEF[0];

    //dudX = K*grad[u] => fraq = (invK*dudX)[0] = (invK*K*grad[u])[0];
    STATE fraq = 0;
    for(int index = 0; index < this->fInvPerm.Cols(); index++)
        fraq += (this->fInvPerm.GetVal(0,index))*dudX[index];

    fraq = fraq - du_exact(0,0);
    REAL diff = fabs(fraq);
	values[3] = diff*diff;
	if(fDim > 1) {
        fraq = 0;
        for(int index = 0; index < this->fInvPerm.Cols(); index++)
            fraq += (this->fInvPerm.GetVal(1,index))*dudX[index];

        fraq = fraq - du_exact(1,0);
		diff = fabs(fraq);
		values[4] = diff*diff;
		if(fDim > 2) {
            fraq = 0;
            for(int index = 0; index < this->fInvPerm.Cols(); index++)
                fraq += (this->fInvPerm.GetVal(1,index))*dudX[index];

            fraq = fraq - du_exact(2,0);
            diff = fabs(fraq);
			values[5] = diff*diff;
		}
	}
	
	TPZManVector<STATE,3> sol(1),dsol(3,0.);
	Solution(u,dudx,axes,1,sol);
	Solution(u,dudx,axes,2,dsol);
	int id;
	//values[1] : eror em norma L2
    diff = fabs(sol[0] - u_exact[0]);
	values[1]  = diff*diff;

	//values[2] : erro em semi norma H1 : |H1| = K*grad[u] - K*grad[u_exact] = K*(grad[u]-grad[u_exact])
	values[2] = 0.;
	TPZVec<REAL> graduDiff(fDim,0);
    for(id=0; id<fDim; id++)  graduDiff[id] += fabs(dsol[id]-du_exact(id,0));
	for(id=0; id<fDim; id++) {
	    STATE diff =0;
        for(int jd = 0; jd <fDim; jd ++){
            diff += fTensorPerm.GetVal(id,jd)*graduDiff[jd];
        }
		values[2]  += abs(diff*diff);
	}
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

void TPZMatPoisson3d::BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump){
    int numbersol = leftu.size();
    for (int is=0; is<numbersol ; is++) {
        jump[is].Resize(1);
        if(bc.Type() == 0){ //DIRICHLET
            STATE f = bc.Val2()(0,0);
            jump[is][0] = leftu[is][0] - f;
        }
        else{
            jump[is].Fill(0.);
        }
    }
}//method

void TPZMatPoisson3d::ContributeEnergy(TPZVec<REAL> &x,
									   TPZVec<FADFADREAL> &sol,
									   TPZVec<FADFADREAL> &dsol,
									   FADFADREAL &U,
									   REAL weight)
{
	int dim = dsol.NElements()/sol.NElements();
	
	//Equa�o de Poisson
	if(sol.NElements() != 1) PZError << "";
	REAL vartocast = 0.;
#ifdef STATE_COMPLEX
    vartocast = fXf.real();
#else
    vartocast = fXf;
#endif
	U+= sol[0] * FADREAL(weight * vartocast);

	STATE avgK = AVGK();
#ifdef STATE_COMPLEX
    vartocast = avgK.real();
#else
    vartocast = avgK;
#endif
	switch(dim)
	{
		case 1:
			U+=vartocast*(dsol[0] * dsol[0])*FADREAL(weight/2.); // U=((du/dx)^2)/2
			
			break;
		case 2:
			U+=vartocast*(dsol[0] * dsol[0] +
				   dsol[1] * dsol[1])*(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2
			/*Buff  = dsol[0] * dsol[0];
             Buff += dsol[1] * dsol[1];
			 U += Buff * FADREAL(weight/2.); // U=((du/dx)^2+(du/dy)^2)/2*/
			break;
		case 3:
			U+=vartocast*(dsol[0] * dsol[0] + dsol[1] * dsol[1] +
				   dsol[2] * dsol[2])*(weight/2.); // U=((du/dx)^2+(du/dy)^2+(du/dz)^2)/2*/
			/*Buff  = dsol[0] * dsol[0];
             Buff += dsol[1] * dsol[1];
             Buff += dsol[2] * dsol[2];
			 U += Buff * FADREAL(weight/2.); //  U=((du/dx)^2+(du/dy)^2+(du/dz)^2)/2*/
			break;
	}
}

void TPZMatPoisson3d::ContributeBCEnergy(TPZVec<REAL> & x,TPZVec<FADFADREAL> & sol, FADFADREAL &U, REAL weight, TPZBndCond &bc) {
    REAL vartocast = 0.;
#ifdef STATE_COMPLEX
    vartocast = bc.Val2()(0,0).real();
#else
    vartocast = bc.Val2()(0,0);
#endif
	FADFADREAL solMinBC = sol[0] - FADREAL(vartocast);
	
	
	switch (bc.Type()) {
		case 0 :	// Dirichlet condition
			// U += 1/2* Big * weight * Integral((u - u0)^2 dOmega)
			U += (solMinBC * solMinBC) * FADREAL(weight * gBigNumber / 2.);
			break;
		case 1 :	// Neumann condition
			// U -= weight * Integral([g].u dOmega)
			U -= sol[0] * FADREAL(vartocast*weight);
			break;
		case 2 :	// condi�o mista
#ifdef STATE_COMPLEX
            vartocast = bc.Val1()(0,0).real();
#else
            vartocast = bc.Val1()(0,0);
#endif
			// U += 1/2 * weight * Integral(<(u-u0), [g].(u-u0)> dOmega)
			U += ( solMinBC * /*scalar*/ FADREAL(vartocast) * /*matrix oprt*/ solMinBC ) * FADREAL(weight / 2.);
			break;
			
	}
}



void TPZMatPoisson3d::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){

    /* This contribute is designed to solve the following  equation:
     *
     *             FIRST_TERM + SECOND_TERM + THIRD_TERM + VOLUME_L-TERM = RIGHT_TERM + VOLUME_R-TERM
     *
     *     FIRST_TERM = CONVECTION TERM : PLEASE DESCRIBE THIS EQUATION HERE!
     *
     *     SECOND_TERM =  -${< K*G[phi_j].n > [[phi_i]] dr} -sym*${< K*G[phi_i].n > [[phi_j]] dr}
     *
     *     IF (PENALTY_CONSTANT == 0.)
     *          THIRD_TERM = 0.
     *
     *     IF (F_PENALTY_TYPE == E_SOLUTION_PENALTY)
     *          THIRD_TERM =  + ${pen*[[phi_j]]*[[phi_i]] dr}
     *
     *     IF (F_PENALTY_TYPE == E_FLUX_PENALTY)
     *          THIRD_TERM = + ${(1/pen)*[[ K*G[phi_j] ]]*[[ phi_i ]] dr}
     *
     *     IF (F_PENALTY_TYPE == BOTH)
     *          THIRD_TERM = + ${(1/pen)*[[ K*G[phi_j] ]]*[[ phi_i ]] dr} + ${pen*[[phi_j]]*[[phi_i]] dr}
     *
     *     RIGHT_TERM = -sym*${<K*G[pd].n> [[phi_j]] dr}
     *
     *     VOLUME_L-TERM && VOLUME_R-TERM: LOOK AT THE VOLUME CONTRIBUTE.
     *
     *     WHERE:
     *
     *    r    : Gamma;
     *  ${.}   : Integral
     *  [[.]]  : Jump operator
     *   <.>   : average operator
     *    .    : dot multiplication (scalar product)
     *  G[.]   : Gradient of
     *   sym   : Symmetry coeficient: -1 or +1
     * pen,pd,n: Penalty term, dirichlet potential and normal vector, respectively.
     *   V_C   : Convection vector
     */
	
	TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
	TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZFMatrix<REAL> &phiR = dataright.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	
	TPZFNMatrix<660> dphiL, dphiR;
	TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
	TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
	
	int &LeftPOrder=dataleft.p;
	int &RightPOrder=dataright.p;
	
	REAL &faceSize=data.HSize;
	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr,id;
	
	//Convection term
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC * fConvDir[id] * normal[id];
	if(ConvNormal > 0.) {
		for(il=0; il<nrowl; il++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
			}
		}
		for(ir=0; ir<nrowr; ir++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(ir+nrowl,jl) -= weight * ConvNormal * phiR(ir) * phiL(jl);
			}
		}
	} else {
		for(ir=0; ir<nrowr; ir++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(ir+nrowl,jr+nrowl) -= weight * ConvNormal * phiR(ir) * phiR(jr);
			}
		}
		for(il=0; il<nrowl; il++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(il,jr+nrowl) += weight * ConvNormal * phiL(il) * phiR(jr);
			}
		}
	}

	// K ?= {{0,0,0}{0,0,0}{0,0,0}}?
	bool isZero = true;
	for(int in = 0; in <fDim; in++) for(int jn =0; jn < fDim; jn++)
	    if(fTensorPerm.GetVal(in,jn) !=0.) isZero = false;
    if(isZero) return;

    TPZFNMatrix<9,STATE> tensorLeftK, tensorRightK;
    tensorLeftK = fTensorPerm;
    tensorRightK = fTensorPerm;

	//KdphiL = K*dphiL ; KdphiR = K*dphiR
    TPZFNMatrix<660,STATE> KdphiL(3,3,0), KdphiR(3,3,0);
	for(int in = 0; in < dphiL.Rows(); in++) for(int jn = 0; jn < dphiL.Cols(); jn++) for(int kn = 0; kn < dphiL.Rows() ; kn++){
	    KdphiL(in,jn) +=  tensorLeftK.GetVal(in,kn)*dphiL(kn,jn);
	    KdphiR(in,jn) +=  tensorRightK.GetVal(in,kn)*dphiR(kn,jn);
	}

	// SECOND_TERM

	// 1) phi_I_left, phi_J_left
	for(il=0; il<nrowl; il++) {
		STATE KdphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			KdphiLinormal += KdphiL(id,il)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE KdphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				KdphiLjnormal += KdphiL(id,jl)*normal[id];
			}
			ek(il,jl) += (STATE)(weight * ( this->fSymmetry * (0.5)*KdphiLinormal*phiL(jl,0)-(0.5)*KdphiLjnormal*phiL(il,0)));
		}
	}
	
	// 2) phi_I_right, phi_J_right
	for(ir=0; ir<nrowr; ir++) {
		STATE KdphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			KdphiRinormal += KdphiR(id,ir)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE KdphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				KdphiRjnormal += KdphiR(id,jr)*normal[id];
			}
			ek(ir+nrowl,jr+nrowl) += (STATE)(weight * (this->fSymmetry * ((-0.5) * KdphiRinormal * phiR(jr)) + (0.5) * KdphiRjnormal * phiR(ir)));
		}
	}
	
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		STATE KdphiLinormal = 0.;
		for(id=0; id<fDim; id++) {
			KdphiLinormal += KdphiL(id,il)*normal[id];
		}
		for(jr=0; jr<nrowr; jr++) {
			STATE KdphiRjnormal = 0.;
			for(id=0; id<fDim; id++) {
				KdphiRjnormal += KdphiR(id,jr)*normal[id];
			}
			ek(il,jr+nrowl) += (STATE)(weight * (STATE)fSymmetry * (STATE)((-0.5) * KdphiLinormal * phiR(jr)) - (STATE)((0.5) * KdphiRjnormal * phiL(il)));
		}
	}
	
	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		STATE KdphiRinormal = 0.;
		for(id=0; id<fDim; id++) {
			KdphiRinormal += KdphiR(id,ir)*normal[id];
		}
		for(jl=0; jl<nrowl; jl++) {
			STATE KdphiLjnormal = 0.;
			for(id=0; id<fDim; id++) {
				KdphiLjnormal += KdphiL(id,jl)*normal[id];
			}
			ek(ir+nrowl,jl) += (STATE)weight * ((STATE)(fSymmetry * (STATE)((0.5) * KdphiRinormal * phiL(jl))  + (STATE)((0.5) * KdphiLjnormal * phiR(ir))));
		}
	}
	
	if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
	
	if (this->fPenaltyConstant == 0.) return;

    // THIRD_TERM

	//penalty = <A p^2>/h
	STATE averageTraceK = AVGK();
	STATE penalty = fPenaltyConstant * (0.5 * averageTraceK*(STATE)(LeftPOrder*LeftPOrder + RightPOrder*RightPOrder)) / faceSize;
	
	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
		
		// 1) left i / left j
		for(il=0; il<nrowl; il++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr; ir++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl; il++) {
			for(jr=0; jr<nrowr; jr++) {
				ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr; ir++) {
			for(jl=0; jl<nrowl; jl++) {
				ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
			}
		}
	}
	
	if (this->fPenaltyType == EFluxPenalty || this->fPenaltyType == EBoth){
		
		REAL NormalFlux_i = 0.;
		REAL NormalFlux_j = 0.;
		
		// 1) left i / left j
		for(il=0; il<nrowl; il++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiL(id,il)*normal[id];
			}
			for(jl=0; jl<nrowl; jl++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiL(id,jl)*normal[id];
				}
				ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j)*averageTraceK;
			}
		}
		
		// 2) right i / right j
		for(ir=0; ir<nrowr; ir++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiR(id,ir)*normal[id];
			}
			for(jr=0; jr<nrowr; jr++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiR(id,jr)*normal[id];
				}      
				ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j)*averageTraceK;
			}
		}
		
		// 3) left i / right j
		for(il=0; il<nrowl; il++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiL(id,il)*normal[id];
			}
			for(jr=0; jr<nrowr; jr++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiR(id,jr)*normal[id];
				}      
				ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j)*averageTraceK;
			}
		}
		
		// 4) right i / left j
		for(ir=0; ir<nrowr; ir++) {
			NormalFlux_i = 0.;
			for(id=0; id<fDim; id++) {
				NormalFlux_i += dphiR(id,ir)*normal[id];
			}
			for(jl=0; jl<nrowl; jl++) {
				NormalFlux_j = 0.;
				for(id=0; id<fDim; id++) {
					NormalFlux_j += dphiL(id,jl)*normal[id];
				}
				ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j)*averageTraceK;
			}
		}
	}
}

/** Termos de penalidade. */
void TPZMatPoisson3d::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, 
                                            REAL weight,
                                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
	
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder= dataleft.p;
	REAL faceSize=data.HSize;

    //KdphiL = K*dphiL ; KdphiR = K*dphiR
    TPZFNMatrix<660,STATE> KdphiL(3,3,0);
    for(int in = 0; in < dphiL.Rows(); in++) for(int jn = 0; jn < dphiL.Cols(); jn++) for(int kn = 0; kn < dphiL.Rows() ; kn++){
        KdphiL(in,jn) +=  fTensorPerm.GetVal(in,kn)*dphiL(kn,jn);
    }
	
	//  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	REAL ConvNormal = 0.;
	for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];
	switch(bc.Type()) {
		case 0: // DIRICHLET
			
			//Diffusion
			for(il=0; il<nrowl; il++) {
				STATE KdphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					KdphiLinormal += KdphiL(id,il)*normal[id];
				}
				ef(il,0) += (STATE)(weight*KdphiLinormal*fSymmetry)*bc.Val2()(0,0);
				for(jl=0; jl<nrowl; jl++) {
					STATE KdphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						KdphiLjnormal += KdphiL(id,jl)*normal[id];
					}
					ek(il,jl) += (STATE)(weight*(fSymmetry * KdphiLinormal * phiL(jl,0) - KdphiLjnormal * phiL(il,0)));
				}
			}
			
			//Convection
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
					}
				}
			} else {
				for(il=0; il<nrowl; il++) {
					ef(il,0) -= (STATE)(weight * ConvNormal * phiL(il)) * bc.Val2()(0,0);
				}
			}
			
			break;
			
		case 1: // Neumann
			for(il=0; il<nrowl; il++) {
				ef(il,0) += (STATE)(weight*phiL(il,0))*bc.Val2()(0,0);
			}
			break;
			
		case 3: // outflow condition
			if(ConvNormal > 0.) {
				for(il=0; il<nrowl; il++) {
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
					}
				}
			}
			else {
				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
			}
			break;
			
		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}
    if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
	
	if (this->fPenaltyConstant == 0.) return;
	
	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){

        STATE averageTraceK = AVGK();
		nrowl = phiL.Rows(); 
		const STATE penalty = fPenaltyConstant * averageTraceK * (STATE)(POrder * POrder) / faceSize; //Ap^2/h
		REAL outflow = 0.;
		for(il=0; il<fDim; il++) outflow += fC * fConvDir[il] * normal[il];
		
		
		switch(bc.Type()) {
			case 0: // DIRICHLET  
				for(il=0; il<nrowl; il++) {
					ef(il,0) += (STATE)(weight * penalty * phiL(il,0)) * bc.Val2()(0,0);
					for(jl=0; jl<nrowl; jl++) {
						ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
					}
				}
				
				break;
			case 1: // Neumann
				if(outflow > 0.)
				{
					for(il=0; il<nrowl; il++)
					{
						for(jl=0; jl<nrowl; jl++)
						{
							ek(il,jl) += weight * outflow * phiL(il,0) * phiL(jl,0);
						}
					}
				}
				//nothing to be done
				break;
			default:
				PZError << "TPZMatPoisson3d::Wrong boundary condition type\n";
				break;
		}
        
	}
	
}

void TPZMatPoisson3d::InterfaceErrors(TPZVec<REAL> &/*x*/,
                                      TPZVec<STATE> &leftu, TPZFMatrix<STATE> &leftdudx, /* TPZFMatrix<REAL> &leftaxes,*/ 
									  TPZVec<STATE> &rightu, TPZFMatrix<STATE> &rightdudx, /* TPZFMatrix<REAL> &rightaxes,*/ 
                                      TPZVec<STATE> &/*flux*/,
									  TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values, 
									  TPZVec<STATE> normal, STATE elsize) {
	// #warning Metodo nao funcional
	TPZManVector<STATE,3> Lsol(1), Ldsol(3,0.), Rsol(1), Rdsol(3,0.);
	
	TPZFMatrix<REAL> fake_axes(fDim,fDim,0.);  
	
	Solution(leftu,leftdudx,fake_axes,1,Lsol);
	Solution(leftu,leftdudx,fake_axes,2,Ldsol);
	
	Solution(rightu,rightdudx,fake_axes,1,Rsol);
	Solution(rightu,rightdudx,fake_axes,2,Rdsol);
	
#ifdef PZDEBUG
	if ( (leftdudx.Rows() != rightdudx.Rows()) || (leftdudx.Rows() != du_exact.Rows()) ){
		PZError << "TPZMatPoisson3d::InterfaceErrors - Left and right matrices should have" 
	    << endl 
	    << "same sizes in internal boundaries." 
	    << endl;
		exit (-1);
	}
#endif
	
	STATE Ldsolnormal = 0., Rdsolnormal = 0., ExactDNormal = 0.;
	for(int id = 0; id < fDim; id++) {
		Ldsolnormal  += Ldsol[id] * normal[id];
		Rdsolnormal  += Rdsol[id] * normal[id];
		ExactDNormal += du_exact(id, 0) * normal[id];
	}
	
	values.Resize(3);
	STATE aux;
	
	//values[1] : eror em norma L2
	
	//Jump aprox. solution - jump of exact solution i.e. zero
	aux = (Lsol[0] - Rsol[0]);
	
	//*= h ^ -gAlfa
	aux *= pow(elsize, (STATE(-1.)) * gAlfa);
    REAL auxnorm = abs(aux);
	values[1] = auxnorm * auxnorm;
	
	//values[2] : erro em semi norma H1
	values[2] = 0.;
	
	for(int id=0; id<fDim; id++) {
		//Normal gradient average <grad V> = 0.5 * (grad_left.n + grad_right.n)
		aux = STATE(0.5) * (Ldsolnormal + Rdsolnormal);
		//<grad V> - <grad exact> = <grad V> - grad exact
		aux = aux - ExactDNormal;
		//*= h ^ gAlfa
		aux *= pow(elsize, gAlfa);
        auxnorm = abs(aux);
		values[2]  += auxnorm * auxnorm;
	}
	//values[0] : erro em norma H1 <=> norma Energia
	values[0]  = values[1]+values[2];
}

REAL TPZMatPoisson3d::ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
	// residual = -fK Laplac(u) + fC * div(fConvDir*u) - (-fXf)
	STATE fXfLoc = fXf;
	if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(X,res);
		fXfLoc = res[0];
	}
	
	STATE laplacU = (STATE)0;
	STATE divBetaU = (STATE)0;
	if(this->Dimension() == 1){
		laplacU = dsol(1,0);
		divBetaU = (STATE)(this->fC * this->fConvDir[0]) * dsol(0,0);
	}
	if(this->Dimension() == 2){
		laplacU = dsol(2,0);
		divBetaU = (STATE)fC * ( (STATE)fConvDir[0] * dsol(0,0) + (STATE)fConvDir[1] * dsol(1,0) );
	} 

	STATE averageTraceK = AVGK();

	REAL result = abs(-averageTraceK * laplacU + divBetaU - (-fXfLoc));
	return (result*result);
}

STATE TPZMatPoisson3d::AVGK(){
    STATE averageTraceK = 0.;
    for(int in =0; in < fTensorPerm.Rows() ; in ++) averageTraceK += fTensorPerm.GetVal(in,in);
    return abs(averageTraceK)/3;
}

void TPZMatPoisson3d::Write(TPZStream &buf, int withclassid) const{
	TPZMaterial::Write(buf, withclassid);
	buf.Write(&fXf, 1);
	buf.Write(&fDim, 1);

	STATE vec1[3], vec2[3], vec3[3];
	for(int in =0 ; in <3 ; in++){
	    vec1[in] = fTensorPerm.Get(0,in);
        vec2[in] = fTensorPerm.Get(1,in);
        vec3[in] = fTensorPerm.Get(2,in);
    }

	buf.Write(vec1, 3);
    buf.Write(vec2, 3);
    buf.Write(vec3, 3);
	buf.Write(&fC, 1);
	buf.Write(fConvDir, 3);
	buf.Write(&fSymmetry, 1);
	buf.Write(&fSD, 1);
	buf.Write(&fPenaltyConstant,1);
	buf.Write(&gAlfa, 1);
}

void TPZMatPoisson3d::Read(TPZStream &buf, void *context){
	TPZMaterial::Read(buf, context);
	buf.Read(&fXf, 1);
	buf.Read(&fDim, 1);

    STATE vec1[3], vec2[3], vec3[3];
    for(int in =0 ; in <3 ; in++){
        vec1[in] = fTensorPerm.Get(0,in);
        vec2[in] = fTensorPerm.Get(1,in);
        vec3[in] = fTensorPerm.Get(2,in);
    }

    buf.Write(vec1, 3);
    buf.Write(vec2, 3);
    buf.Write(vec3, 3);
	buf.Read(&fC, 1);
	buf.Read(fConvDir, 3);
	buf.Read(&fSymmetry, 1);
	buf.Read(&fSD, 1);
	buf.Read(&fPenaltyConstant,1);
	buf.Read(&gAlfa, 1);
}

int TPZMatPoisson3d::ClassId() const{
    return Hash("TPZMatPoisson3d") ^ TPZMaterial::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZMatPoisson3d>;
#endif
