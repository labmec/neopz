/**
 * @file
 * @brief Contains implementations of the TPZMatDarcyMHM methods.
 */

#include "DarcyMHM.h"
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include "pzmaterialdata.h"
#include <math.h>
#include "pzlog.h"

#include <cmath>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.poisson3d"));
#endif


using namespace std;

TPZMatDarcyMHM::TPZMatDarcyMHM(int nummat, int dim) : TPZDiscontinuousGalerkin(nummat), fXf(0.), fDim(dim), fMultiplier(1.) {
	fK = 1.;
}

TPZMatDarcyMHM::TPZMatDarcyMHM():TPZDiscontinuousGalerkin(), fXf(0.), fDim(1), fMultiplier(1.){
	fK = 1.;
}

TPZMatDarcyMHM::TPZMatDarcyMHM(const TPZMatDarcyMHM &copy):TPZDiscontinuousGalerkin(copy){
	this->operator =(copy);
}

TPZMatDarcyMHM & TPZMatDarcyMHM::operator=(const TPZMatDarcyMHM &copy){
	TPZDiscontinuousGalerkin::operator = (copy);
	fXf  = copy.fXf;
	fDim = copy.fDim;
	fK   = copy.fK;
    fMultiplier = copy.fMultiplier;
	return *this;
}


TPZMatDarcyMHM::~TPZMatDarcyMHM() {
}

int TPZMatDarcyMHM::NStateVariables() {
	return 1;
}

void TPZMatDarcyMHM::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Laplace operator multiplier fK "<< fK << endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}

void TPZMatDarcyMHM::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef) {
	
    //return;
	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &dphi = data.dphix;
	TPZVec<REAL>  &x = data.x;
	TPZFMatrix<REAL> &axes = data.axes;
	TPZFMatrix<REAL> &jacinv = data.jacinv;
	int phr = phi.Rows();

	STATE fXfLoc = 1;//fXf;
	
//	if(fForcingFunction)
//    {            // phi(in, 0) = phi_in
//		TPZManVector<STATE,1> res(1);
//		TPZFMatrix<STATE> dres(Dimension(),1);
//		fForcingFunction->Execute(x,res,dres);       // dphi(i,j) = dphi_j/dxi
//		fXfLoc = res[0];
//	}
//	REAL delx = 0.;

	for( int in = 0; in < phr; in++ )
    {
		int kd;
		ef(in, 0) += weight * fXfLoc * ( (STATE)phi(in,0)  );
		for( int jn = 0; jn < phr; jn++ )
        {
			for(kd=0; kd<fDim; kd++)
            {
				ek(in,jn) += (STATE)weight * fK * (STATE)( dphi(kd,in) * dphi(kd,jn) );
			}
		}
	}
    
    
}


void TPZMatDarcyMHM::ContributeBC(TPZMaterialData &data,REAL weight,
								   TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc) {
	
	//return;
	TPZFMatrix<REAL>  &phi = data.phi;
	TPZFMatrix<REAL> &axes = data.axes;
	int phr = phi.Rows();
	short in,jn;
	STATE v2[1];
	v2[0] = bc.Val2()(0,0);
	
//	if(bc.HasForcingFunction())            // phi(in, 0) = phi_in
//    {
//		TPZManVector<STATE> res(1);
//		bc.ForcingFunction()->Execute(data.x,res);       // dphi(i,j) = dphi_j/dxi
//		v2[0] = res[0];
//	}

	switch (bc.Type())
    {
		case 0 :			// Dirichlet condition
			for(in = 0 ; in < phr; in++)
            {
				ef(in,0) += (STATE)(gBigNumber* phi(in,0) * weight) * v2[0];
				for (jn = 0 ; jn < phr; jn++)
                {
					ek(in,jn) += gBigNumber * phi(in,0) * phi(jn,0) * weight;
				}
			}
			break;
		case 1 :			// Neumann condition
//			for(in = 0 ; in < phi.Rows(); in++)
//            {
//				ef(in,0) += v2[0] * (STATE)(phi(in,0) * weight);
//			}
			break;
		case 2 :		// mixed condition
//			for(in = 0 ; in < phi.Rows(); in++)
//            {
//				ef(in, 0) += v2[0] * (STATE)(phi(in, 0) * weight);
//				for (jn = 0 ; jn < phi.Rows(); jn++)
//                {
//					ek(in,jn) += bc.Val1()(0,0) * (STATE)(phi(in,0) * phi(jn,0) * weight);     // peso de contorno => integral de contorno
//				}
//			}
			break;
		case 3: // outflow condition
//			int id, il, jl;
//			REAL normal[3];
//			if (fDim == 1) PZError << __PRETTY_FUNCTION__ << " - ERROR! The normal vector is not available for 1D TPZInterpolatedElement\n";
//			if (fDim == 2)
//            {
//				normal[0] = axes(0,1);
//				normal[1] = axes(1,1);
//			}
//			if (fDim == 3)
//            {
//				normal[0] = axes(0,2);
//				normal[1] = axes(1,2);
//				normal[2] = axes(2,2);
//			}
//			REAL ConvNormal = 0.;    
////			for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];  
//			if(ConvNormal > 0.)
//            {
//				for(il=0; il<phr; il++)
//                {
//					for(jl=0; jl<phr; jl++)
//                    {
//						ek(il,jl) += weight * ConvNormal * phi(il)*phi(jl);
//					}
//				}
//			}
//			else{
//				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
//			}
			break;
	}
	}

/** Returns the variable index associated with the name */
int TPZMatDarcyMHM::VariableIndex(const std::string &name){
	if(!strcmp("Solution",name.c_str()))        return  1;
//	if(!strcmp("Derivative",name.c_str()))      return  2;
//	if(!strcmp("KDuDx",name.c_str()))           return  3;
//	if(!strcmp("KDuDy",name.c_str()))           return  4;
//	if(!strcmp("KDuDz",name.c_str()))           return  5;
//	if(!strcmp("NormKDu",name.c_str()))         return  6;
//	if(!strcmp("MinusKGradU",name.c_str()))     return  7;
//	if(!strcmp("POrder",name.c_str()))          return  8;
//	if(!strcmp("Laplac",name.c_str()))          return  9;
//	if(!strcmp("Stress",name.c_str()))          return  10;    
//	if(!strcmp("Flux",name.c_str()))            return  10;
//	if(!strcmp("Pressure",name.c_str()))        return  11;
//	
//	if(!strcmp("ExactPressure",name.c_str()))   return  12;
//	if(!strcmp("ExactFlux",name.c_str()))       return  13;
//	if(!strcmp("Divergence",name.c_str()))      return  14;
//	if(!strcmp("ExactDiv",name.c_str()))        return  15;
//	
//	if(!strcmp("PressureOmega1",name.c_str()))  return  16;
//	if(!strcmp("PressureOmega2",name.c_str()))  return  17;
//	if(!strcmp("FluxOmega1",name.c_str()))      return  18;
//    
//    if(!strcmp("GradFluxX",name.c_str()))       return  19;
//    if(!strcmp("GradFluxY",name.c_str()))       return  20;
//     if(!strcmp("FluxL2",name.c_str()))            return  21;//Only To calculate l2 error
	return TPZMaterial::VariableIndex(name);
}

int TPZMatDarcyMHM::NSolutionVariables(int var){
	if(var == 1) return 1;
//	if(var == 2) return fDim;//arrumar o fluxo de hdiv para ser fdim tbem enquanto isso faco isso
//	if ((var == 3) || (var == 4) || (var == 5) || (var == 6)) return 1;
//	if (var == 7) return fDim;
//	if (var == 8) return 1;
//	if (var == 9) return 1;
//	if (var==10) return fDim;
//	if (var==11) return 1;
//	
//	if (var==12) return 1;
//	if (var==13) return fDim;
//	if (var==14) return 1;
//	if (var==15) return 1;
//	//teste de acoplamento
//	if (var==16) return 1;
//	if (var==17) return 1;
//	if (var==18) return 3;
//    if (var==19) return 3;
//    if (var==20) return 3;
//    if (var==21) return fDim;
	
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMatDarcyMHM::Solution(TPZMaterialData &data, int var, TPZVec<STATE> &Solout){
	
	TPZVec<STATE> pressure(1);
	TPZVec<REAL> pto(3);
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
            Solout[0]=data.sol[0][data.sol[0].NElements()-1];
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
            this->Solution(data.sol[0], data.dsol[0], data.axes, var, Solout);
            break;
    }
#endif
}

#include "pzaxestools.h"
void TPZMatDarcyMHM::Solution(TPZVec<STATE> &Sol,TPZFMatrix<STATE> &DSol,TPZFMatrix<REAL> &axes,int var,TPZVec<STATE> &Solout){
	
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


void TPZMatDarcyMHM::Errors(TPZVec<REAL> &x,TPZVec<STATE> &u,
							 TPZFMatrix<STATE> &dudx, TPZFMatrix<REAL> &axes, TPZVec<STATE> &/*flux*/,
							 TPZVec<STATE> &u_exact,TPZFMatrix<STATE> &du_exact,TPZVec<REAL> &values){
    return;
//	values.Resize(NEvalErrors());
//	TPZManVector<STATE> dudxEF(1,0.), dudyEF(1,0.),dudzEF(1,0.);
//	this->Solution(u,dudx,axes,VariableIndex("KDuDx"), dudxEF);
//    STATE fraq = dudxEF[0]/fK;
//    fraq = fraq - du_exact(0,0);
//    REAL diff = fabs(fraq);
//	values[3] = diff*diff;
//	if(fDim > 1) {
//		this->Solution(u,dudx,axes, this->VariableIndex("KDuDy"), dudyEF);
//        fraq = dudyEF[0]/fK;
//        fraq = fraq - du_exact(1,0);
//		diff = fabs(fraq);
//		values[4] = diff*diff;
//		if(fDim > 2) {
//			this->Solution(u,dudx,axes, this->VariableIndex("KDuDz"), dudzEF);
//			fraq = dudzEF[0]/fK;
//            fraq = fraq - du_exact(2,0);
//            diff = fabs(fraq);
//			values[5] = diff*diff;
//		}
//	}
//	
//	TPZManVector<STATE> sol(1),dsol(3,0.);
//	Solution(u,dudx,axes,1,sol);
//	Solution(u,dudx,axes,2,dsol);
//	int id;
//	//values[1] : eror em norma L2
//    diff = fabs(sol[0] - u_exact[0]);
//	values[1]  = diff*diff;
//	//values[2] : erro em semi norma H1
//	values[2] = 0.;
//	for(id=0; id<fDim; id++) {
//        diff = fabs(dsol[id] - du_exact(id,0));
//		values[2]  += abs(fK)*diff*diff;
//	}
//	//values[0] : erro em norma H1 <=> norma Energia
//	values[0]  = values[1]+values[2];
}




void TPZMatDarcyMHM::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
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
	
	
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr;
    
	// 3) phi_I_left, phi_J_right
	for(il=0; il<nrowl; il++) {
		for(jr=0; jr<nrowr; jr++) {
			ek(il,jr+nrowl) += weight * fMultiplier * (phiL(il) * phiR(jr));
		}
	}
	
//	// 4) phi_I_right, phi_J_left
	for(ir=0; ir<nrowr; ir++) {
		for(jl=0; jl<nrowl; jl++) {
			ek(ir+nrowl,jl) += weight * fMultiplier * (phiR(ir) * phiL(jl));
		}
	}
		
	
}

void TPZMatDarcyMHM::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft,
                                            REAL weight,
                                            TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    int teste=0;
    return;
	
	TPZFMatrix<REAL> &dphiL = dataleft.dphix;
	TPZFMatrix<REAL> &phiL = dataleft.phi;
	TPZManVector<REAL,3> &normal = data.normal;
	int POrder= dataleft.p;
	REAL faceSize=data.HSize;
	
	//  cout << "Material Id " << bc.Id() << " normal " << normal << "\n";
	int il,jl,nrowl,id;
	nrowl = phiL.Rows();
	REAL ConvNormal = 0.;
//	for(id=0; id<fDim; id++) ConvNormal += fC*fConvDir[id]*normal[id];
	switch(bc.Type()) {
		case 0: // DIRICHLET
			
			//Diffusion
			for(il=0; il<nrowl; il++) {
				REAL dphiLinormal = 0.;
				for(id=0; id<fDim; id++) {
					dphiLinormal += dphiL(id,il)*normal[id];
				}
//				ef(il,0) += (STATE)(weight*dphiLinormal*fSymmetry)*fK*bc.Val2()(0,0);
				for(jl=0; jl<nrowl; jl++) {
					REAL dphiLjnormal = 0.;
					for(id=0; id<fDim; id++) {
						dphiLjnormal += dphiL(id,jl)*normal[id];
					}
//					ek(il,jl) += (STATE)(weight*(fSymmetry * dphiLinormal * phiL(jl,0) - dphiLjnormal * phiL(il,0)))*fK;
				}
			}
            break;
			
		case 1: // Neumann
//			for(il=0; il<nrowl; il++)
//            {
//				ef(il,0) += (STATE)(weight*phiL(il,0))*bc.Val2()(0,0);
//			}
			break;
			
		case 3: // outflow condition
//			if(ConvNormal > 0.)
//            {
//				for(il=0; il<nrowl; il++)
//                {
//					for(jl=0; jl<nrowl; jl++)
//                    {
//						ek(il,jl) += weight * ConvNormal * phiL(il)*phiL(jl);
//					}
//				}
//			}
//			else {
//				if (ConvNormal < 0.) std::cout << "Boundary condition error: inflow detected in outflow boundary condition: ConvNormal = " << ConvNormal << "\n";
//			}
			break;
			
		default:
			PZError << __PRETTY_FUNCTION__ << " - Wrong boundary condition type\n";
			break;
	}

	
}


void TPZMatDarcyMHM::Write(TPZStream &buf, int withclassid){
	TPZDiscontinuousGalerkin::Write(buf, withclassid);
	buf.Write(&fXf, 1);
	buf.Write(&fDim, 1);
	buf.Write(&fK, 1);
}

void TPZMatDarcyMHM::Read(TPZStream &buf, void *context){
	TPZDiscontinuousGalerkin::Read(buf, context);
	buf.Read(&fXf, 1);
	buf.Read(&fDim, 1);
	buf.Read(&fK, 1);
}

