/**
 * @file
 * @brief Contains implementations of the tpzdifureac methods
 */

#include "tpzdifureac.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"
#include "pzaxestools.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.tpzdifureac"));// esto no sé para qué es
#endif




using namespace std;

TPZdifureac::TPZdifureac(int nummat, int dim) : TPZDiscontinuousGalerkin(nummat), fXf(1.), fDim(dim){ //constructor
  fK = 0.; /////////////////////
  falpha = 0;
  fPenaltyConstant = 1000.;
  this->SetNonSymmetric();
  this->SetNoPenalty();
}

TPZdifureac::TPZdifureac():TPZDiscontinuousGalerkin(), fXf(0.), fDim(1){
	fK = 0.;
	falpha = 0;
	fPenaltyConstant = 1000.;
	this->SetNonSymmetric();
	this->SetNoPenalty();
}

TPZdifureac::TPZdifureac(const TPZdifureac &copy):TPZDiscontinuousGalerkin(copy){
	this->operator =(copy);//const. de copia
}


TPZdifureac & TPZdifureac::operator=(const TPZdifureac &copy){
	TPZDiscontinuousGalerkin::operator = (copy);
	fXf  = copy.fXf;
	fDim = copy.fDim;
	fK   = copy.fK;
	falpha = copy.falpha;
	fSymmetry = copy.fSymmetry;
	fPenaltyConstant = copy.fPenaltyConstant;
	this->fPenaltyType = copy.fPenaltyType;
	return *this;
}

void TPZdifureac::SetParameters(STATE diff, STATE f, STATE alf) {
  fK = diff;
  fXf = f;
  falpha = alf;
}

void TPZdifureac::GetParameters(STATE &diff, STATE &f, STATE &alf) {
    diff = fK;
    f = fXf;
    alf = falpha;
}

TPZdifureac::~TPZdifureac()
{
//destrutor
}


void TPZdifureac::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Laplace operator multiplier fK "<< fK << endl;
	out << "Forcing vector fXf " << fXf << endl;
	out << "Reaction term: "<<falpha<< endl;
	out << "Penalty constant " << fPenaltyConstant << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out); 
	out << "\n";
}

void TPZdifureac::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	
       
    TPZFMatrix<REAL>  &phi = data.phi;
    TPZFMatrix<REAL> &dphi = data.dphix;
    TPZVec<REAL>  &x = data.x;
    //    TPZFMatrix<REAL> &axes = data.axes;
    //    TPZFMatrix<REAL> &jacinv = data.jacinv;
    int phr = phi.Rows();
    
    STATE fXfLoc = fXf;
    
    
    if(fForcingFunction) {            // phi(in, 0) = phi_in
        TPZManVector<STATE,1> res(1);	//se tiver termo fonte, calcula no veto
        TPZFMatrix<STATE> dres(Dimension(),1);
        fForcingFunction->Execute(x,res);       // dphi(i,j) = dphi_j/dxi
        fXfLoc = res[0];
    }
    
    //Equacao de Poisson
    for( int in = 0; in < phr; in++ ) { // e independente da dimensão
        int kd;
        ef(in, 0) +=  (STATE)weight * fXfLoc * (STATE)phi(in,0);
        
	 for( int jn = 0; jn < phr; jn++ ) {
	  ek(in,jn) += (STATE)weight * (phi(in,0)*phi(jn,0))*falpha; 
	   for(kd=0; kd<fDim; kd++) {
                ek(in,jn) += (STATE)weight * (
                                              +fK * (STATE)( dphi(kd,in) * dphi(kd,jn)) ); //
            }
        }
    }
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() )
	{
	  cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
    }
    
}


void TPZdifureac::ContributeBC(TPZMaterialData &data,REAL weight,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
}

/** Returns the variable index associated with the name */
int TPZdifureac::VariableIndex(const std::string &name){
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
    
	if(!strcmp("ExactPressure",name.c_str()))   return  12;
	if(!strcmp("ExactFlux",name.c_str()))       return  13;
	if(!strcmp("Divergence",name.c_str()))      return  14;
	if(!strcmp("ExactDiv",name.c_str()))        return  15;
    
	if(!strcmp("PressureOmega1",name.c_str()))  return  16;
	if(!strcmp("PressureOmega2",name.c_str()))  return  17;
	if(!strcmp("FluxOmega1",name.c_str()))      return  18;
    
    if(!strcmp("GradFluxX",name.c_str()))       return  19;
    if(!strcmp("GradFluxY",name.c_str()))       return  20;
    if(!strcmp("FluxL2",name.c_str()))            return  21;//Only To calculate l2 error
    if(!strcmp("ExactSolution",name.c_str()))        return  22;
    if(!strcmp("ExactDerivative",name.c_str()))        return  23;    
	return TPZMaterial::VariableIndex(name);
}

int TPZdifureac::NSolutionVariables(int var){
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
    if (var==22) return 1;
    if (var==23) return fDim;     
    
    
	return TPZMaterial::NSolutionVariables(var);
}

void TPZdifureac::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
    
	TPZVec<STATE> pressure(1);
	TPZVec<REAL> pto(3);
	TPZVec<STATE> SolExact(1);	
	TPZFMatrix<STATE> flux(3,1);
    
    int numbersol = data.sol.size();
    if (numbersol != 1) {
        DebugStop();
    }
    
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
		case 22:
		{
		  fForcingFunctionExact->Execute(data.x,SolExact,flux);
		  Solout[0] = SolExact[0];
		}
			break;	
		case 23:
		{
		  fForcingFunctionExact->Execute(data.x,SolExact,flux);
		  Solout[0] = flux[0];
		  Solout[1] = flux[1];
		}
			break;				
		case 10:
			if (data.numberdualfunctions) {
                
				Solout[0]=data.sol[0][0];
				Solout[1]=data.sol[0][1];
                
			}
			else {
				this->Solution(data.sol[0], data.dsol[0], data.axes, 2, Solout);
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
		case 13:
            fForcingFunctionExact->Execute(data.x,pressure,flux);
            
            Solout[0]=flux(0,0);
            Solout[1]=flux(1,0);
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
        default:
            if (data.sol[0].size() == 4) {
                data.sol[0][0] = data.sol[0][2];
            }
            
            this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
            break;
    }
#endif
}

void TPZdifureac::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
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

void TPZdifureac::ErrorsHdiv(TPZMaterialData &data,TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){
	
	TPZVec<STATE> sol(1),dsol(fDim),div(1);
	if(data.numberdualfunctions) Solution(data,11,sol);//pressao
	Solution(data,21,dsol);//fluxo
	Solution(data,14,div);//divergente
    
#ifdef LOG4CXX
    {
		std::stringstream sout;
		sout<< "\n";
		sout << " Pto  " << data.x << std::endl;
		sout<< " pressao exata " <<u_exact <<std::endl;
		sout<< " pressao aprox " <<sol <<std::endl;
		sout<< " ---- "<<std::endl;
		sout<< " fluxo exato " <<du_exact<<std::endl;
		sout<< " fluxo aprox " <<dsol<<std::endl;
		sout<< " ---- "<<std::endl;
		if(du_exact.Rows()>fDim) sout<< " div exato " <<du_exact(2,0)<<std::endl;
		sout<< " div aprox " <<div<<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
    }
#endif
    
    
	//values[0] : pressure error using L2 norm
	if(data.numberdualfunctions){
		REAL diffP = abs(u_exact[0]-sol[0]);
		values[0]  = diffP*diffP;
	}	
	//values[1] : flux error using L2 norm
	for(int id=0; id<fDim; id++) {
        REAL diffFlux = abs(dsol[id] - du_exact(id,0));
		values[1]  += abs(fK)*diffFlux*diffFlux;
	}
    
    if(du_exact.Rows()>fDim){
        //values[2] : divergence using L2 norm
        REAL diffDiv = abs(div[0] - du_exact(2,0));
        values[2]=diffDiv*diffDiv;
        //values[3] : Hdiv norm => values[1]+values[2];
        values[3]= values[1]+values[2];
    }
    
#ifdef LOG4CXX
    {
		std::stringstream sout;
		sout << " Erro pressao  " << values[0]<< std::endl;
		sout<< "Erro fluxo  " <<values[1]<<std::endl;
		sout<< " Erro div " <<values[2] <<std::endl;
		LOGPZ_DEBUG(logger,sout.str())
    }
#endif
	
    
}

void TPZdifureac::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
							 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes,
							 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values) {
	

	
	TPZFNMatrix<9,STATE> dsoldx;
	TPZAxesTools<STATE>::Axes2XYZ(dudx, dsoldx, axes);
  
	values.Resize(3);
	///L2 norm
	values[1] = (u[0] - u_exact[0])*(u[0] - u_exact[0]);
	
	///semi norma de H1
	values[2] = 0.;
	for(int i = 0; i < this->fDim; i++){
		values[2] += (dsoldx(i,0) - du_exact(i,0))*(dsoldx(i,0) - du_exact(i,0));	
	}
	
	///H1 norm
	values[0] = values[1]+values[2]; //donde saca la raíz ?
	
}

void TPZdifureac::BCInterfaceJump(TPZVec<REAL> &x, TPZSolVec &leftu,TPZBndCond &bc,TPZSolVec & jump){
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


void TPZdifureac::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                          REAL weight,
                                          TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef){
    
    TPZFMatrix<REAL> &dphiLdAxes = dataleft.dphix;
    TPZFMatrix<REAL> &dphiRdAxes = dataright.dphix;
    TPZFMatrix<REAL> &phiL = dataleft.phi;
    TPZFMatrix<REAL> &phiR = dataright.phi;
    TPZManVector<REAL,3> &normal = data.normal;
    
    REAL sign = 1.0; 
    
    normal[0] *= sign;
    normal[1] *= sign;
    normal[2] *= sign;    
    
    TPZFNMatrix<660> dphiL, dphiR;
    TPZAxesTools<REAL>::Axes2XYZ(dphiLdAxes, dphiL, dataleft.axes);
    TPZAxesTools<REAL>::Axes2XYZ(dphiRdAxes, dphiR, dataright.axes);
    
//    int &LeftPOrder=dataleft.p;
//    int &RightPOrder=dataright.p;
    
    REAL &faceSize=data.HSize;
    
    
    int nrowl = phiL.Rows();
    int nrowr = phiR.Rows();
    int il,jl,ir,jr,id;
   
    REAL penaljumpuv = this->fsigma/pow(faceSize,this->fbeta);
    
    //diffusion term
    STATE leftK, rightK;
    leftK  = this->fK;
    rightK = this->fK;
    
    // 1) phi_I_left, phi_J_left M11
    for(il=0; il<nrowl; il++) {
        REAL dphiLinormal = 0.;//gradphi*normal
        for(id=0; id<fDim; id++) {
            dphiLinormal += dphiL(id,il)*normal[id];
        }
        for(jl=0; jl<nrowl; jl++) {
            REAL dphiLjnormal = 0.;
            for(id=0; id<fDim; id++) {
                dphiLjnormal += dphiL(id,jl)*normal[id];
            }
            ek(il,jl) += (STATE)(weight * ( (this->fSymmetry * (0.5)*dphiLinormal*phiL(jl,0)-(0.5)*dphiLjnormal*phiL(il,0)) * leftK +
            penaljumpuv*phiL(jl,0)*phiL(il,0) ) );
        }
    }
    // 2) phi_I_right, phi_J_right M22
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
            ek(ir+nrowl,jr+nrowl) += (STATE)(weight * (this->fSymmetry * ((-0.5) * dphiRinormal * phiR(jr) )*rightK + (0.5) * dphiRjnormal * phiR(ir) * rightK 
            + penaljumpuv*phiR(jr)*phiR(ir) ) );
        }
    }
    
    // 3) phi_I_left, phi_J_right M12
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
            ek(il,jr+nrowl) += (STATE)(weight *( this->fSymmetry * (-0.5) * dphiLinormal * phiR(jr) * leftK 
            -((0.5) * dphiRjnormal * phiL(il))*rightK - penaljumpuv*phiR(jr)*phiL(il)  ) );
        }
    }
    
    // 4) phi_I_right, phi_J_left M21
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
            ek(ir+nrowl,jl) += (STATE) ( weight * ( this->fSymmetry *(0.5) *dphiRinormal*phiL(jl)*rightK + (0.5)*dphiLjnormal*phiR(ir)*leftK
           -penaljumpuv*phiL(jl)*phiR(ir)) );
        }
    }
    
    if (this->IsSymetric()){
        if ( !ek.VerifySymmetry() )
	{
	  cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
	}
    }
    
//     if (this->fPenaltyConstant == 0.) return;
//     
//     leftK  = this->fK;
//     rightK = this->fK;
//     
//     
//     
//     //penalty = <A p^2>/h
//     REAL penalty = fPenaltyConstant * (0.5 * (abs(leftK)*LeftPOrder*LeftPOrder + abs(rightK)*RightPOrder*RightPOrder)) / faceSize;
//     
//     if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){
//         
//         // 1) left i / left j
//         for(il=0; il<nrowl; il++) {
//             for(jl=0; jl<nrowl; jl++) {
//                 ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
//             }
//         }
//         
//         // 2) right i / right j
//         for(ir=0; ir<nrowr; ir++) {
//             for(jr=0; jr<nrowr; jr++) {
//                 ek(ir+nrowl,jr+nrowl) += weight * penalty * phiR(ir,0) * phiR(jr,0);
//             }
//         }
//         
//         // 3) left i / right j
//         for(il=0; il<nrowl; il++) {
//             for(jr=0; jr<nrowr; jr++) {
//                 ek(il,jr+nrowl) += -1.0 * weight * penalty * phiR(jr,0) * phiL(il,0);
//             }
//         }
//         
//         // 4) right i / left j
//         for(ir=0; ir<nrowr; ir++) {
//             for(jl=0; jl<nrowl; jl++) {
//                 ek(ir+nrowl,jl) += -1.0 * weight *  penalty * phiL(jl,0) * phiR(ir,0);
//             }
//         }
//         
//     }
//     
//     if (this->fPenaltyType == EFluxPenalty || this->fPenaltyType == EBoth){
//         
//         REAL NormalFlux_i = 0.;
//         REAL NormalFlux_j = 0.;
//         
//         // 1) left i / left j
//         for(il=0; il<nrowl; il++) {
//             NormalFlux_i = 0.;
//             for(id=0; id<fDim; id++) {
//                 NormalFlux_i += dphiL(id,il)*normal[id];
//             }
//             for(jl=0; jl<nrowl; jl++) {
//                 NormalFlux_j = 0.;
//                 for(id=0; id<fDim; id++) {
//                     NormalFlux_j += dphiL(id,jl)*normal[id];
//                 }
//                 ek(il,jl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
//             }
//         }
//         
//         // 2) right i / right j
//         for(ir=0; ir<nrowr; ir++) {
//             NormalFlux_i = 0.;
//             for(id=0; id<fDim; id++) {
//                 NormalFlux_i += dphiR(id,ir)*normal[id];
//             }
//             for(jr=0; jr<nrowr; jr++) {
//                 NormalFlux_j = 0.;
//                 for(id=0; id<fDim; id++) {
//                     NormalFlux_j += dphiR(id,jr)*normal[id];
//                 }
//                 ek(ir+nrowl,jr+nrowl) += (STATE)(weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
//             }
//         }
//         
//         // 3) left i / right j
//         for(il=0; il<nrowl; il++) {
//             NormalFlux_i = 0.;
//             for(id=0; id<fDim; id++) {
//                 NormalFlux_i += dphiL(id,il)*normal[id];
//             }
//             for(jr=0; jr<nrowr; jr++) {
//                 NormalFlux_j = 0.;
//                 for(id=0; id<fDim; id++) {
//                     NormalFlux_j += dphiR(id,jr)*normal[id];
//                 }
//                 ek(il,jr+nrowl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * rightK;
//             }
//         }
//         
//         // 4) right i / left j
//         for(ir=0; ir<nrowr; ir++) {
//             NormalFlux_i = 0.;
//             for(id=0; id<fDim; id++) {
//                 NormalFlux_i += dphiR(id,ir)*normal[id];
//             }
//             for(jl=0; jl<nrowl; jl++) {
//                 NormalFlux_j = 0.;
//                 for(id=0; id<fDim; id++) {
//                     NormalFlux_j += dphiL(id,jl)*normal[id];
//                 }
//                 ek(ir+nrowl,jl) += (STATE)((-1.) * weight * ((1.)/penalty) * NormalFlux_i * NormalFlux_j) * leftK;
//             }
//         }
//         
//     }
//     
}

/** Termos de penalidade. */
void TPZdifureac::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                            REAL weight,
                                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	REAL faceSize=data.HSize;
	REAL penaljumpuv = 2*this->fsigma/pow(faceSize,this->fbeta);
	
	REAL sign = 1.0; 
    
	normal[0] *= sign;
	normal[1] *= sign;
	normal[2] *= sign;    

	
	//  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	switch(bc.Type()) {
		case 0: // DIRICHLET
            {
			//Diffusion
			for(il=0; il<nrowl; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL(id,il)*normal[id];
				}
				ef(il,0) += (STATE)(weight * (dphiLinormal*fSymmetry*fK+ penaljumpuv*phiL(il,0))*bc.Val2()(0,0));
  
				for(jl=0; jl<nrowl; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL(id,jl)*normal[id];
					}
					
					ek(il,jl) += (STATE)( weight*(
					  (fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0)) * fK
					  + penaljumpuv*phiL(jl,0)*phiL(il,0)));
				}
			}
	    }
	    break;
		case 1: // Neumann
		{
            for(il=0; il<nrowl; il++) {
				ef(il,0) += (STATE)(weight*phiL(il,0))*bc.Val2()(0,0);
			}
		}
			break;
            
		default:
		{
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
		}
			break;
	}
    if (this->IsSymetric()){
		if ( !ek.VerifySymmetry() ) 
		{
		  cout << __PRETTY_FUNCTION__ << "\nMATRIZ NAO SIMETRICA" << endl;
		}
    }
    
		
//  if (this->fenaltyConstant == 0.) return;
//     
// 	if (this->fPenaltyType == ESolutionPenalty || this->fPenaltyType == EBoth){ //o qué faz isto?
// 		nrowl = phiL.Rows();
// 		const REAL penalty = fPenaltyConstant * abs(fK) * POrder * POrder / faceSize; //Ap^2/h
// 		REAL outflow = 0.;
//         
// 		switch(bc.Type()) {
// 			case 0: // DIRICHLET
// 				for(il=0; il<nrowl; il++) {
// 					ef(il,0) += (STATE)(weight * penalty * phiL(il,0)) * bc.Val2()(0,0);
// 					for(jl=0; jl<nrowl; jl++) {
// 						ek(il,jl) += weight * penalty * phiL(il,0) * phiL(jl,0);
// 					}
// 				}
//                 
// 				break;
// 			case 1: // Neumann
// 				if(outflow > 0.)
// 				{
// 					for(il=0; il<nrowl; il++)
// 					{
// 						for(jl=0; jl<nrowl; jl++)
// 						{
// 							ek(il,jl) += weight * outflow * phiL(il,0) * phiL(jl,0);
// 						}
// 					}
// 				}
// 				//nothing to be done
// 				break;
// 			default:
// 				PZError << "TPZdifureac::Wrong boundary condition type\n";
// 				break;
// 		}
//         
   
}

REAL TPZdifureac::ComputeSquareResidual(TPZVec<REAL>& X, TPZVec<STATE> &sol, TPZFMatrix<STATE> &dsol){
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

void TPZdifureac::Write(TPZStream &buf, int withclassid) const{
	TPZDiscontinuousGalerkin::Write(buf, withclassid);
	buf.Write(&fXf, 1);
	buf.Write(&fDim, 1);
	buf.Write(&fK, 1);
	buf.Write(&fSymmetry, 1);
	buf.Write(&fPenaltyConstant,1);
}

void TPZdifureac::Read(TPZStream &buf, void *context){
	TPZDiscontinuousGalerkin::Read(buf, context);
	buf.Read(&fXf, 1);
	buf.Read(&fDim, 1);
	buf.Read(&fK, 1);
	buf.Read(&fSymmetry, 1);
	buf.Read(&fPenaltyConstant,1);
}

int TPZdifureac::ClassId() const{
    return Hash("TPZdifureac") ^ TPZDiscontinuousGalerkin::ClassId() << 1;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZdifureac>; //no sé qué deba ser aqui
#endif
