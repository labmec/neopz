/**
 * @file
 * @brief Contains implementations of the TPZMatLaplacian methods.
 */

#include "TPZMatLaplacian.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"
#include "pzextractval.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.TPZMatLaplacian"));
#endif


using namespace std;

TPZMatLaplacian::TPZMatLaplacian(int nummat, int dim) :
TPZRegisterClassId(&TPZMatLaplacian::ClassId), TPZMaterial(nummat), fXf(0.), fDim(dim), fTensorK(dim,dim,0.), fInvK(dim,dim,0.)
{
	fK = 1.;
    for (int i=0; i<dim; i++) {
        fTensorK(i,i) = 1.;
        fInvK(i,i) = 1.;
    }
	fPenaltyConstant = 1000.;
	this->SetNonSymmetric();
	this->SetNoPenalty();
}

TPZMatLaplacian::TPZMatLaplacian()
: TPZRegisterClassId(&TPZMatLaplacian::ClassId), TPZMaterial(), fXf(0.), fDim(1), fTensorK(1,1,1.), fInvK(1,1,1.){
	fK = 1.;
	fPenaltyConstant = 1000.;
	this->SetNonSymmetric();
	this->SetNoPenalty();
}

TPZMatLaplacian::TPZMatLaplacian(const TPZMatLaplacian &copy)
: TPZRegisterClassId(&TPZMatLaplacian::ClassId), TPZMaterial(copy)
{
	this->operator =(copy);
}

TPZMatLaplacian & TPZMatLaplacian::operator=(const TPZMatLaplacian &copy){
	TPZMaterial::operator = (copy);
	fXf  = copy.fXf;
	fDim = copy.fDim;
	fK   = copy.fK;
    fTensorK = copy.fTensorK;
    fInvK = copy.fInvK;
	fSymmetry = copy.fSymmetry;
	fPenaltyConstant = copy.fPenaltyConstant;
	this->fPenaltyType = copy.fPenaltyType;
    this->fPermeabilityFunction = copy.fPermeabilityFunction;
	return *this;
}

void TPZMatLaplacian::SetParameters(STATE diff, STATE f) {
	fK = diff;
    fTensorK.Zero();
    fInvK.Zero();
    for (int i=0; i<fDim; i++) {
        fTensorK(i,i) = diff;
        fInvK(i,i) = 1./diff;
    }
    fXf = f;
}

void TPZMatLaplacian::SetPermeabilityTensor(const TPZFNMatrix<9,STATE> K, const TPZFNMatrix<9,STATE> invK){

    if(K.Rows() != 3 || invK.Rows() != 3)
    {
        if(K.Rows() != 2 || invK.Rows() != 2) {
            std::cout << "ERROR: Insert a 3x3 or a 2x2 permeability tensor";
            DebugStop();
        }
        for(int iind =0; iind < 2 ; iind++) for(int jind = 0; jind <2 ; jind++){
                fTensorK(iind,jind) = K.GetVal(iind,jind);
                fInvK(iind,jind) = invK.GetVal(iind,jind);
            }
    }
    else {
        fTensorK = K;
        fInvK = invK;
    }
}

void TPZMatLaplacian::GetPermeability(TPZFNMatrix<9,STATE> &K){
    K = fTensorK;
}

void TPZMatLaplacian::GetInvPermeability(TPZFNMatrix<9,STATE> &invK){
    invK = fInvK;
}

void TPZMatLaplacian::GetPermeabilities(TPZVec<REAL> &x, TPZFNMatrix<9,STATE> &PermTensor, TPZFNMatrix<9,STATE> &InvPermTensor)
{
    this->GetPermeability(PermTensor);
    this->GetInvPermeability(InvPermTensor);
}

TPZMatLaplacian::~TPZMatLaplacian() {
}

int TPZMatLaplacian::NStateVariables() const {
	return 1;
}

void TPZMatLaplacian::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Laplace operator multiplier fK "<< fK << endl;
	out << "Forcing vector fXf " << fXf << endl;
	out << "Penalty constant " << fPenaltyConstant << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMatLaplacian::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    STATE KPerm = fK;
    if (fPermeabilityFunction) {
        TPZFNMatrix<9,STATE> perm, invperm;
        TPZManVector<STATE,3> func;
        TPZFNMatrix<18,STATE> dfunc(6,3,0.);
        fPermeabilityFunction->Execute(x, func, dfunc);
        KPerm = dfunc(0,0);
    }

    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for( int jn = 0; jn < phr; jn++ ) {
            //ek(in,jn) += (STATE)weight*((STATE)(phi(in,0)*phi(jn,0)));
            for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight*(KPerm*(STATE)(dphi(kd,in)*dphi(kd,jn)));
            }
        }
    }
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry(1.e-10) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
    
}

void TPZMatLaplacian::Contribute(TPZMaterialData &data,REAL weight, TPZFMatrix<STATE> &ef)
{
    
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);
        //TPZFMatrix<STATE> dres(Dimension(),1);
        //fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
        fForcingFunction->Execute(x,res);
        fXfLoc = res[0];
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) {
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        for(kd=0; kd<fDim; kd++) {
            ef(in,0) -= (STATE)weight*(fK*(STATE)(dphi(kd,in)*data.dsol[0](kd,0)));
        }
    }
}

void TPZMatLaplacian::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {

    ContributeBC(datavec[0],weight,ek,ef,bc);
}

void TPZMatLaplacian::ContributeBC(TPZMaterialData &data,REAL weight,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
	TPZFMatrix<REAL>  &phi = data.phi;
    //	TPZFMatrix<REAL> &axes = data.axes;
	int phr = phi.Rows();
	short in,jn;
	STATE v2[1];
	v2[0] = bc.Val2()(0,0);
    
	if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in                          // JORGE 2013 01 26
		TPZManVector<STATE> res(1);
        TPZFNMatrix<3,STATE> dres(3,1);
		bc.ForcingFunction()->Execute(data.x,res,dres);       // dphi(i,j) = dphi_j/dxi
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
	}
    
	if (this->IsSymetric()) {//only 1.e-3 because of bignumbers.
		if ( !ek.VerifySymmetry( 1.e-3 ) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
}

/** Returns the variable index associated with the name */
int TPZMatLaplacian::VariableIndex(const std::string &name){
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
    
	if(!strcmp("ExactSolution",name.c_str()))   return  12;
	if(!strcmp("ExactFlux",name.c_str()))       return  13;
	if(!strcmp("Divergence",name.c_str()))      return  14;
	if(!strcmp("ExactDiv",name.c_str()))        return  15;
    
	if(!strcmp("PressureOmega1",name.c_str()))  return  16;
	if(!strcmp("PressureOmega2",name.c_str()))  return  17;
	if(!strcmp("FluxOmega1",name.c_str()))      return  18;
    
    if(!strcmp("GradFluxX",name.c_str()))       return  19;
    if(!strcmp("GradFluxY",name.c_str()))       return  20;
    if(!strcmp("FluxL2",name.c_str()))          return  21;//Only To calculate l2 error
    if(!strcmp("Permeability",name.c_str()))    return  22; // output the permeability

    if(!strcmp("ExactFluxShiftedOrigin",name.c_str()))  return 23; //Shift the coordinates at origin by 10^-10

	return TPZMaterial::VariableIndex(name);
}

int TPZMatLaplacian::NSolutionVariables(int var){
	if(var == 1) return 1;
	if(var == 2) return 3;//arrumar o fluxo de hdiv para ser fdim tbem enquanto isso faco isso
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
    if (var==22) return 1; // number of permeabilities

    if (var == 23) return fDim;
    
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatLaplacian::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
	TPZVec<STATE> pressure(1);
	TPZVec<REAL> pto(3);
	TPZFMatrix<STATE> flux(3,1);
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
    STATE perm = fK;
    if(fPermeabilityFunction)
    {
        TPZManVector<STATE,3> f;
        TPZFNMatrix<18,STATE> df(6,3);
        fPermeabilityFunction->Execute(data.x, f, df);
        perm = df(0,0);
    }

    // Solution EArcTan returns NAN for (x,y) = (0,0). Replacing data.x by inf solves this problem,
    STATE infinitesimal = 0.0000000001;
    TPZManVector<STATE,3> inf ={infinitesimal,infinitesimal,infinitesimal};

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
                
			}
			else {
                TPZFNMatrix<3,STATE> dsolxy(3,0);
                TPZAxesTools<STATE>::Axes2XYZ(data.dsol[0], dsolxy, data.axes);
                for (int i=0; i<fDim; i++) {
                    Solout[i] = -perm*dsolxy(i,0);
                }
			}
            
			break;
            
        case 21:
            Solout[0]=data.sol[0][0];
            Solout[1]=data.sol[0][1];
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
            fForcingFunctionExact->Execute(data.x,pressure,flux);
            
            Solout[0]=pressure[0];
			break;
		case 13: //ExactFlux
		    fForcingFunctionExact->Execute(data.x, pressure, flux);

            Solout[0]=-flux(0,0);
            Solout[1]=-flux(1,0);
            break;

        case 23: //ExactFluxShiftedOrigin
            if(data.x[0] == 0. && data.x[1] == 0.) {
                fForcingFunctionExact->Execute(inf, pressure, flux);
            } else {
                fForcingFunctionExact->Execute(data.x, pressure, flux);
            }

            if (std::isnan(flux(0, 0))) {
                std::cout << "Flux X is NAN at: " << data.x[0] << "," << data.x[1] << std::endl;
            }
            if (std::isnan(flux(1, 0))) {
                std::cout << "Flux Y is NAN at: " << data.x[0] << "," << data.x[1] << std::endl;
            }

            Solout[0]=-flux(0,0);
            Solout[1]=-flux(1,0);
            break;
            
        case 14:
        {
			if (data.numberdualfunctions){
				Solout[0]=data.sol[0][data.sol[0].NElements()-1];
			}else{
				Solout[0]=data.dsol[0](0,0)+data.dsol[0](1,1);
			}
            
        }
            break;
            
        case 15:
        {
            fForcingFunctionExact->Execute(data.x,pressure,flux);
            Solout[0]=flux(2,0);
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
        case 22:
            // output the permeability
            Solout[0] = perm;
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

void TPZMatLaplacian::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
#ifndef STATE_COMPLEX
	Solout.Resize( this->NSolutionVariables( var ) );
	
	if(var == 1){
		Solout[0] = Sol[0];//function
		return;
	}
	if(var == 2) {
		int id;
        TPZFNMatrix<50,STATE> dsoldx;
        TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		for(id=0 ; id<3; id++) {
			Solout[id] = dsoldx(id,0);//derivate
		}
		return;
	}//var == 2
	if (var == 3){ //KDuDx
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(0,0) * this->fK;
		return;
	}//var ==3
	if (var == 4){ //KDuDy
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(1,0) * this->fK;
		return;
	}//var == 4 
	if (var == 5){ //KDuDz
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		Solout[0] = dsoldx(2,0) * this->fK;
		return;
	}//var == 5
	if (var == 6){ //NormKDu
		int id;
		REAL val = 0.;
		for(id=0 ; id<fDim; id++){
			val += (DSol(id,0) * this->fK) * (DSol(id,0) * this->fK);
		}
		Solout[0] = sqrt(val);
		return;
	}//var == 6
	if (var == 7){ //MinusKGradU
		int id;
		//REAL val = 0.;
		TPZFNMatrix<9,STATE> dsoldx;
		TPZAxesTools<STATE>::Axes2XYZ(DSol, dsoldx, axes);
		for(id=0 ; id<fDim; id++) {
			Solout[id] = -1. * this->fK * dsoldx(id,0);
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

void TPZMatLaplacian::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
							 TPZFMatrix<STATE> &dudxaxes, TPZFMatrix<REAL> &axes, 
							 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	
    values.Resize(3);
    values.Fill(0.0);

    TPZFNMatrix<9,STATE> perm(2*fDim,fDim);
    TPZManVector<STATE,3> val(fDim);
    if (fPermeabilityFunction) {
        fPermeabilityFunction->Execute(x, val, perm);
    }
    else
    {
        for (int i=0; i<fDim; i++) {
            for (int j=0; j<fDim; j++)
            {
                perm(i,j) = this->fTensorK(i,j);
                perm(fDim+i,j) = this->fInvK(i,j);
            }
        }
    }
    TPZFNMatrix<3,STATE> dudx(3,1,0.);
    TPZAxesTools<STATE>::Axes2XYZ(dudxaxes, dudx, axes);
    
	///L2 norm
	values[1] = TPZExtractVal::val((u[0] - u_exact[0])*(u[0] - u_exact[0]));
	
	///semi norma de H1
	values[2] = 0.;
	for(int i = 0; i < du_exact.Rows(); i++){
		values[2] += TPZExtractVal::val( (dudx(i,0) - du_exact(i,0))*(dudx(i,0) - du_exact(i,0)));
	}
    // Energy Norm
    values[0] = 0.;
    for (int i=0; i<fDim; i++) {
        for (int j=0; j<fDim; j++) {
                values[0] += TPZExtractVal::val((dudx(i,0) - du_exact(i,0))*perm(i,j)*(dudx(j,0) - du_exact(j,0)));
        }
    }
	///H1 norm

}

void TPZMatLaplacian::BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump){
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


void TPZMatLaplacian::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
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
    
    //diffusion term
    STATE leftK, rightK;
    leftK  = this->fK;
    rightK = this->fK;
    
    // 1) phi_I_left, phi_J_left
    for(il=0; il<nrowl; il++) {
        REAL dphiLinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiLinormal += dphiL(id,il)*normal[id];
        }
        for(jl=0; jl<nrowl; jl++) {
            REAL dphiLjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiLjnormal += dphiL(id,jl)*normal[id];
            }
            ek(il,jl) += (STATE)(weight * ( this->fSymmetry * (0.5)*dphiLinormal*phiL(jl,0)-(0.5)*dphiLjnormal*phiL(il,0))) * leftK;
        }
    }
    
    // 2) phi_I_right, phi_J_right
    for(ir=0; ir<nrowr; ir++) {
        REAL dphiRinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiRinormal += dphiR(id,ir)*normal[id];
        }
        for(jr=0; jr<nrowr; jr++) {
            REAL dphiRjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiRjnormal += dphiR(id,jr)*normal[id];
            }
            ek(ir+nrowl,jr+nrowl) += (STATE)(weight * (this->fSymmetry * ((-0.5) * dphiRinormal * phiR(jr) ) + (0.5) * dphiRjnormal * phiR(ir))) * rightK;
        }
    }
    
    // 3) phi_I_left, phi_J_right
    for(il=0; il<nrowl; il++) {
        REAL dphiLinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiLinormal += dphiL(id,il)*normal[id];
        }
        for(jr=0; jr<nrowr; jr++) {
            REAL dphiRjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiRjnormal += dphiR(id,jr)*normal[id];
            }
            ek(il,jr+nrowl) += (STATE)weight * ((STATE)fSymmetry * ((STATE)((-0.5) * dphiLinormal * phiR(jr)) * leftK ) - (STATE)((0.5) * dphiRjnormal * phiL(il))* rightK );
        }
    }
    
    // 4) phi_I_right, phi_J_left
    for(ir=0; ir<nrowr; ir++) {
        REAL dphiRinormal = 0.;
        for(id=0; id<fDim; id++) {
            dphiRinormal += dphiR(id,ir)*normal[id];
        }
        for(jl=0; jl<nrowl; jl++) {
            REAL dphiLjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiLjnormal += dphiL(id,jl)*normal[id];
            }
            ek(ir+nrowl,jl) += (STATE)weight * (
                                                (STATE)(fSymmetry * (0.5) * dphiRinormal * phiL(jl)) * rightK + (STATE)((0.5) * dphiLjnormal * phiR(ir)) * leftK
                                                );
        }
    }
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry(1.e-10) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
    
    if (this->fPenaltyConstant == 0.) return;
    
    leftK  = this->fK;
    rightK = this->fK;
    
    
    
    //penalty = <A p^2>/h
    REAL penalty = fPenaltyConstant * (0.5 * (abs(leftK)*LeftPOrder*LeftPOrder + abs(rightK)*RightPOrder*RightPOrder)) / faceSize;
    
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
                ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
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
                ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
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
                ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
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
                ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
            }
        }
        
    }
    
}

/** Termos de penalidade. */
void TPZMatLaplacian::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                            REAL weight,
                                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder= dataleft.p;
	REAL faceSize=data.HSize;
    
    STATE v2[1];
    v2[0] = bc.Val2()(0,0);
    
    if(bc.HasForcingFunction()) {            // phi(in, 0) = phi_in
        TPZManVector<STATE> res(1);
        bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
        v2[0] = res[0];
    }
    
	//  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	switch(bc.Type()) {
		case 0: // DIRICHLET
            
			//Diffusion
			for(il=0; il<nrowl; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL(id,il)*normal[id];
				}
				ef(il,0) += (STATE)(weight*dphiLinormal*fSymmetry)*fK*v2[0];
				for(jl=0; jl<nrowl; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL(id,jl)*normal[id];
					}
					ek(il,jl) += (STATE)(weight*(fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0)))*fK;
				}
			}
            
		case 1: // Neumann
			for(il=0; il<nrowl; il++) {
				ef(il,0) += (STATE)(weight*phiL(il,0))*v2[0];
			}
			break;
            
		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}
    if (this->IsSymetric()){
		if ( !ek.VerifySymmetry(1.e-10) ) cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
    }
    
	if (this->fPenaltyConstant == 0.) return;
    
	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
		nrowl = phiL.Rows();
		const REAL penalty = fPenaltyConstant * abs(fK) * POrder * POrder / faceSize; //Ap^2/h
		REAL outflow = 0.;
        
		switch(bc.Type()) {
			case 0: // DIRICHLET
				for(il=0; il<nrowl; il++) {
					ef(il,0) += (STATE)(weight*penalty*phiL(il,0))*v2[0];
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
				PZError << "TPZMatLaplacian::Wrong boundary condition type\n";
				break;
		}
        
	}
    
}


void TPZMatLaplacian::FillDataRequirements(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsNeighborSol = false;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
    data.fNeedsSol = true;
}

void TPZMatLaplacian::FillDataRequirementsInterface(TPZMaterialData &data)
{
    data.SetAllRequirements(false);
    data.fNeedsNeighborSol = false;
    data.fNeedsNeighborCenter = false;
    data.fNeedsNormal = true;
    data.fNeedsHSize = true;
}

REAL TPZMatLaplacian::ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
	// residual = -fK Laplac(u) - (fXf)
	STATE fXfLoc = fXf;
	if(fForcingFunction) {
		TPZManVector<STATE> res(1);
		fForcingFunction->Execute(X,res);
		fXfLoc = res[0];
	}
    
	STATE laplacU = (STATE)0;
	if(this->Dimension() == 1){
		laplacU = dsol(1,0);
	}
	if(this->Dimension() == 2){
		laplacU = dsol(2,0);
	}
    
	REAL result = abs(-this->fK * laplacU - (fXfLoc));
	return (result*result);
}

void TPZMatLaplacian::Write(TPZStream &buf, int withclassid) const {
	TPZMaterial::Write(buf, withclassid);
	buf.Write(&fXf, 1);
	buf.Write(&fDim, 1);
	buf.Write(&fK, 1);
	buf.Write(&fSymmetry, 1);
	buf.Write(&fPenaltyConstant,1);
}

void TPZMatLaplacian::Read(TPZStream &buf, void *context) {
	TPZMaterial::Read(buf, context);
	buf.Read(&fXf, 1);
	buf.Read(&fDim, 1);
	buf.Read(&fK, 1);
	buf.Read(&fSymmetry, 1);
	buf.Read(&fPenaltyConstant,1);
}

int TPZMatLaplacian::ClassId() const{
    return Hash("TPZMatLaplacian") ^ TPZMaterial::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZMatLaplacian>;
#endif
