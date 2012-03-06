//$Id: pzelasAXImat.cpp,v 1.10 2011-03-29 10:56:17 phil Exp $
/**
 * \file
 * @brief Contains implementations of the TPZElasticityAxiMaterial methods.
 */
#include "pzelasAXImat.h" 
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <cmath>
#include "pzaxestools.h"
#include "pzlog.h"

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.material.axisymetric"));
static LoggerPtr logdata(Logger::getLogger("pz.material.axisymetric.data"));
#endif

#include <fstream>
using namespace std;

TPZElasticityAxiMaterial::TPZElasticityAxiMaterial() : TPZDiscontinuousGalerkin(0), fIntegral(0.), fAlpha(1.e-5), f_AxisR(3,0.), f_AxisZ(3,0.),f_Origin(3,0.), fTemperatureFunction(0) {
	f_AxisZ[1] = 1.;
	f_AxisR[0] = 1.;
    fDelTemperature = 0.;
	fE	= -1.;  // Young modulus
	fnu	= -1.;   // poisson coefficient
	ff[0]	= 0.; // X component of the body force
	ff[1]	= 0.; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = -1.;  //G = E/2(1-nu);
	fEover21PlusNu = -1.;//E/(1-nu)
	
	f_c = 0.;
	f_phi = 0.;
    fSymmetric = 1.;
    fPenalty = 1.;
}

TPZElasticityAxiMaterial::TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy) : TPZDiscontinuousGalerkin(num), fIntegral(0.), fAlpha(1.e-5), fDelTemperature(0.), f_AxisR(3,0.), f_AxisZ(3,0.),f_Origin(3,0.), fTemperatureFunction(0)
{
	
    f_AxisZ[1] = 1.;
    f_AxisR[0] = 1.;
	
	fE	= E;  // Young modulus
	fnu	= nu;   // poisson coefficient
	ff[0]	= fx; // X component of the body force
	ff[1]	= fy; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = E/(1-fnu*fnu);  //G = E/2(1-nu);
	fEover21PlusNu = E/(2.*(1+fnu));//E/(1-nu)
	f_c = 0.;
	f_phi = 0.;
    fSymmetric = 1.;
    fPenalty = 1.;
}

//--------------------------------------------------------------------------------------------------------------------------------------
TPZElasticityAxiMaterial::TPZElasticityAxiMaterial(int num, REAL E, REAL nu, REAL fx, REAL fy, REAL coefTheta, REAL coefAlpha) : 
TPZDiscontinuousGalerkin(num), fIntegral(0.), fAlpha(1.e-5), fDelTemperature(0.), f_AxisR(3,0.), f_AxisZ(3,0.),f_Origin(3,0.),
fTemperatureFunction(0)
{
	
	f_AxisZ[1] = 1.;
	f_AxisR[0] = 1.;
	fE	= E;  // Young modulus
	fnu	= nu;   // poisson coefficient
	ff[0]	= fx; // X component of the body force
	ff[1]	= fy; // Y component of the body force
	ff[2] = 0.; // Z component of the body force - not used for this class
	fEover1MinNu2 = E/(1-fnu*fnu);  //G = E/2(1-nu);
	fEover21PlusNu = E/(2.*(1+fnu));//E/(1-nu)
	f_c = 0.;
	f_phi = 0.;
	fSymmetric = coefTheta;
	fPenalty = coefAlpha;
}

/*
 REAL fIntegral;
 REAL f_phi;
 REAL f_c;
 REAL fE;
 REAL fnu;
 REAL fAlpha;
 REAL ff[3];
 REAL fDelTemperature;
 REAL fEover21PlusNu;
 REAL fEover1MinNu2;
 TPZManVector<REAL> f_AxisR;
 TPZManVector<REAL> f_AxisZ;
 TPZManVector<REAL> f_Origin;
 REAL fSymmetric;
 REAL fPenalty;
 */


TPZElasticityAxiMaterial::TPZElasticityAxiMaterial(const TPZElasticityAxiMaterial &copy) : 
TPZDiscontinuousGalerkin(copy), fIntegral(copy.fIntegral),f_phi(copy.f_phi),f_c(copy.f_c), fE(copy.fE),
fnu(copy.fnu), fAlpha(copy.fAlpha), fDelTemperature(copy.fDelTemperature), fEover21PlusNu(copy.fEover21PlusNu),
fEover1MinNu2(copy.fEover1MinNu2),f_AxisR(copy.f_AxisR),f_AxisZ(copy.f_AxisZ),
f_Origin(copy.f_Origin),fSymmetric(copy.fSymmetric),fPenalty(copy.fPenalty),fTemperatureFunction(copy.fTemperatureFunction)
{
	ff[0] = copy.ff[0];
	ff[1] = copy.ff[1];
	ff[2] = copy.ff[2];
}

//--------------------------------------------------------------------------------------------------------------------------------------

void TPZElasticityAxiMaterial::SetOrigin(TPZManVector<REAL> &Orig, TPZManVector<REAL> &AxisZ, TPZManVector<REAL> &AxisR)
{
	if(Orig.size() == 3 && AxisZ.size() == 3 && AxisR.size() == 3)
	{
		TPZFNMatrix<6> Vecs(3,2,0.), VecsOrt(3,2,0.), JacVecsOrth(2,2,0.);
		for(int i = 0; i < 3; i++)
		{
			Vecs.PutVal(i,0,AxisR[i]);
			Vecs.PutVal(i,1,AxisZ[i]);
		}
		Vecs.GramSchmidt(VecsOrt,JacVecsOrth);
		
		for(int j = 0; j < 3; j++)
		{
			AxisR[j] = VecsOrt.GetVal(j,0);
			AxisZ[j] = VecsOrt.GetVal(j,1);
		}
		f_Origin = Orig;
		f_AxisZ = AxisZ;
		f_AxisR = AxisR;
	}
	else
	{
		cout << "Invalid Origin and/or Axis vector on TPZElasticityAxiMaterial()!\n";
		DebugStop();
	}
}

REAL TPZElasticityAxiMaterial::ComputeR(TPZVec<REAL> &x)
{
	return (x[0] - f_Origin[0])*f_AxisR[0] + (x[1] - f_Origin[1])*f_AxisR[1] + (x[2] - f_Origin[2])*f_AxisR[2];
}

TPZManVector<REAL> TPZElasticityAxiMaterial::GetAxisR()
{
	return f_AxisR;
}

TPZManVector<REAL> TPZElasticityAxiMaterial::GetAxisZ()
{
	return f_AxisZ;
}

TPZManVector<REAL> TPZElasticityAxiMaterial::GetOrigin()
{
	return f_Origin;
}

TPZElasticityAxiMaterial::~TPZElasticityAxiMaterial() {
}

int TPZElasticityAxiMaterial::NStateVariables() {
	return 2;
}

void TPZElasticityAxiMaterial::Print(std::ostream &out) {
	TPZMaterial::Print(out);
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\tE   = " << fE   << endl;
	out << "\tnu   = " << fnu   << endl;
	out << "\tF   = " << ff[0] << ' ' << ff[1]   << endl;
}

void TPZElasticityAxiMaterial::Contribute(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef)
{
	TPZFMatrix &dphi = data.dphix;
	TPZFMatrix &phi  = data.phi;
	TPZFMatrix &axes = data.axes;
	
	int phc,phr,dphc,dphr,efr,efc,ekr,ekc;
	phc = phi.Cols();
	phr = phi.Rows();
	dphc = dphi.Cols();
	dphr = dphi.Rows();
	efr = ef.Rows();
	efc = ef.Cols();
	ekr = ek.Rows();
	ekc = ek.Cols();
	if(phc != 1 || dphr != 2 || phr != dphc || ekr != phr*2 || ekc != phr*2 || efr != phr*2 || efc != 1)
	{
		PZError << "\nTPZElasticityMaterial.contr, inconsistent input data : \n" <<
		"phi.Cols() = " << phi.Cols() << " dphi.Cols() = " << dphi.Cols() <<
		" phi.Rows = " << phi.Rows() << " dphi.Rows = " <<
		dphi.Rows() << "\nek.Rows() = " << ek.Rows() << " ek.Cols() = "
	    << ek.Cols() <<
		"\nef.Rows() = " << ef.Rows() << " ef.Cols() = "
	    << ef.Cols() << "\n";
		return;
	}
	if(fForcingFunction)
	{
		TPZManVector<REAL> res(3);
		fForcingFunction->Execute(data.x,res);
		ff[0] = res[0];
		ff[1] = res[1];
		ff[2] = res[2];
	}
	
	// R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!
	REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];
    REAL Z = (data.x[0] - f_Origin[0])*f_AxisZ[0] + (data.x[1] - f_Origin[1])*f_AxisZ[1] + (data.x[2] - f_Origin[2])*f_AxisZ[2];
	
	int s = (R > 0)? 1:-1;
	R = fabs(R);
	if(R < 1.e-10) R = 1.e-10;
	
    REAL DelTemp(fDelTemperature);
	if (fTemperatureFunction) {
        TPZManVector<REAL,2> RZ(2);
        RZ[0] = R;
        RZ[1] = Z;
        fTemperatureFunction(RZ,DelTemp);
    }
	/**
	 * Plain strain materials values
	 */
	// REAL nu1 = 1. - fnu;//(1-nu)
	//  REAL nu2 = (1.-2.*fnu)/2.;
	//  REAL F = fE/((1.+fnu)*(1.-2.*fnu));
	
	TPZFNMatrix<4> dphiRZi(2,1), dphiRZj(2,1);
	
	REAL axis0DOTr = 0., axis1DOTr = 0., axis0DOTz = 0., axis1DOTz = 0.;
	for(int pos = 0; pos < 3; pos++)
	{
		axis0DOTr += axes.GetVal(0,pos) * f_AxisR[pos] * s;
		axis1DOTr += axes.GetVal(1,pos) * f_AxisR[pos] * s;
		axis0DOTz += axes.GetVal(0,pos) * f_AxisZ[pos];
		axis1DOTz += axes.GetVal(1,pos) * f_AxisZ[pos];
	}
	
	REAL R2PI = 2. * M_PI * R;
    
    REAL lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
    REAL mi =  fE/(2.*(1. + fnu));
    
    REAL epsT = DelTemp*fAlpha;
    REAL TensorThermico = (3.*lambda+2.*mi)*epsT;
	
	//criado para resolver o problema de Girkmann	
	fIntegral += R2PI*weight;
	
	for( int in = 0; in < phr; in++ )
	{
		
		//dphi_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
		dphiRZi.PutVal(0,0, dphi.GetVal(0,in)*axis0DOTr + dphi.GetVal(1,in)*axis1DOTr );
		
		//dphi_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
		dphiRZi.PutVal(1,0, dphi.GetVal(0,in)*axis0DOTz + dphi.GetVal(1,in)*axis1DOTz );
		
		ef(2*in, 0)   += weight * R2PI * (ff[0] * phi(in,0) + TensorThermico*(phi(in,0)/R + dphiRZi(0,0))); // direcao x
		ef(2*in+1, 0) += weight * R2PI * (ff[1] * phi(in,0) + TensorThermico*dphiRZi(1,0)); // direcao y
		
		for( int jn = 0; jn < phr; jn++ )
		{
			
			//dphi_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
			dphiRZj.PutVal(0,0, dphi.GetVal(0,jn)*axis0DOTr + dphi.GetVal(1,jn)*axis1DOTr );
			
			//dphi_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
			dphiRZj.PutVal(1,0, dphi.GetVal(0,jn)*axis0DOTz + dphi.GetVal(1,jn)*axis1DOTz );
			
			REAL term00 = dphiRZi(0,0) * (lambda + 2.*mi) * dphiRZj(0,0) +
			dphiRZi(1,0) * mi * dphiRZj(1,0) +
			phi(in,0) * lambda/R * dphiRZj(0,0) +
			dphiRZi(0,0) * lambda/R * phi(jn,0) +
			phi(in,0) * (lambda+2.*mi)/(R*R) * phi(jn,0);
			
			REAL term01 = dphiRZi(1,0) * mi * dphiRZj(0,0) +
			dphiRZi(0,0) * lambda * dphiRZj(1,0) +
			phi(in,0) * lambda/R * dphiRZj(1,0);
			
			REAL term10 = dphiRZi(1,0) * lambda * dphiRZj(0,0) +
			dphiRZi(0,0) * mi * dphiRZj(1,0) +
			dphiRZi(1,0) * lambda/R * phi(jn,0);
			
			REAL term11 = dphiRZi(0,0) * mi * dphiRZj(0,0) +
			dphiRZi(1,0) * (lambda + 2.*mi) * dphiRZj(1,0);
			
			ek(2*in,2*jn)     += weight * R2PI * term00;
			ek(2*in,2*jn+1)   += weight * R2PI * term01;
			ek(2*in+1,2*jn)   += weight * R2PI * term10;
			ek(2*in+1,2*jn+1) += weight * R2PI * term11;
		}
	}
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
}

void TPZElasticityAxiMaterial::ContributeBC(TPZMaterialData &data,REAL weight,TPZFMatrix &ek,TPZFMatrix &ef,TPZBndCond &bc)
{
	TPZFMatrix &phi = data.phi;
	
	const REAL BIGNUMBER  = TPZMaterial::gBigNumber;
	
	int phr = phi.Rows();
	short in,jn;
	REAL v2[2];
	v2[0] = bc.Val2()(0,0);
	v2[1] = bc.Val2()(1,0);
	
	// R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!
	REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];
	
	int s = (R > 0) ? 1:-1;
	R = fabs(R);
	double R2PI = 2. * M_PI * R;
	
	static REAL accum1 = 0., accum2 = 0.;
	
	switch (bc.Type())
	{
		case 0 :// Dirichlet condition
		{
			for(in = 0 ; in < phr; in++)
			{
				ef(2*in,0)   += BIGNUMBER * v2[0]*  // x displacement
                phi(in,0) * R2PI * weight;        // forced v2 displacement
				
				ef(2*in+1,0) += BIGNUMBER * v2[1]*// x displacement
                phi(in,0) * R2PI * weight;        // forced v2 displacement
				
				for (jn = 0 ; jn < phi.Rows(); jn++)
				{
					ek(2*in,2*jn)     += BIGNUMBER * phi(in,0) * phi(jn,0) * R2PI * weight;
					ek(2*in+1,2*jn+1) += BIGNUMBER * phi(in,0) * phi(jn,0) * R2PI * weight;
				}
			}
		}
			break;
			
		case 1 :// Neumann condition
		{
			for(in = 0 ; in < phi.Rows(); in++)
			{           // componentes da tracao na direcao de v2
				ef(2*in,0)   += v2[0] * phi(in,0) * R2PI * weight;   // tracao em x  (ou pressao)
				ef(2*in+1,0) += v2[1] * phi(in,0) * R2PI * weight; // tracao em y (ou pressao) , nula se n� h
			}      // ou deslocamento nulo  v2 = 0
			accum1 += v2[1] * R2PI *weight;
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Accumulated force 1 " << accum1;
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 2 :// condicao mista
		{
			for(in = 0 ; in < phi.Rows(); in++)
			{
				ef(2*in, 0)   += v2[0] * phi(in, 0) * R2PI * weight;   // Neumann , Sigmaij
				ef(2*in+1, 0) += v2[1] * phi(in, 0) * R2PI * weight; // Neumann
				
				for (jn = 0 ; jn < phi.Rows(); jn++)
				{
					ek(2*in,2*jn)     += bc.Val1()(0,0) * phi(in,0) * phi(jn,0) * R2PI * weight; // peso de contorno => integral de contorno
					ek(2*in+1,2*jn)   += bc.Val1()(1,0) * phi(in,0) * phi(jn,0) * R2PI * weight;
					ek(2*in+1,2*jn+1) += bc.Val1()(1,1) * phi(in,0) * phi(jn,0) * R2PI * weight;
					ek(2*in,2*jn+1)   += bc.Val1()(0,1) * phi(in,0) * phi(jn,0) * R2PI * weight;
				}
			}   // este caso pode reproduzir o caso 0 quando o deslocamento
		}
			break;
			
		case 3 :// Neumann condition - Normal rotacionada 90 graus no sentido horário em relacao ao vetor axes (1D)
		{
			REAL Nxy[2],Nrz[2];
			Nxy[0] =  data.axes(0,1);
			Nxy[1] = -data.axes(0,0);
			Nrz[0] = s*(Nxy[0]*f_AxisR[0] + Nxy[1]*f_AxisR[1]);
			Nrz[1] = Nxy[0]*f_AxisZ[0] + Nxy[1]*f_AxisZ[1];
			
			//           #ifdef LOG4CXX
			//           {
			//               std::stringstream sout;
			//               sout << "CoordX: " << data.x << " R= " << R << " Nrz = " << Nrz[0] << "," << Nrz[1] << endl;
			//               LOGPZ_DEBUG(logger,sout.str());
			//           }
			//           #endif
			
			for(in = 0 ; in < phi.Rows(); in++)
			{           // componentes da tracao normal ao contorno
				ef(2*in,0)   += v2[0] * Nrz[0] * phi(in,0) * R2PI * weight;   // tracao em x  (ou pressao)
				ef(2*in+1,0) += v2[0] * Nrz[1] * phi(in,0) * R2PI * weight; // tracao em y (ou pressao) , nula se n� h
			}      // ou deslocamento nulo  v2 = 0
			accum2 += v2[0] * Nrz[1] * R2PI *weight;
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Accumulated force 2 " << accum2;
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
	}      // �nulo introduzindo o BIGNUMBER pelos valores da condicao
}         // 1 Val1 : a leitura �00 01 10 11

//---------------------------- --------------------------------------
void TPZElasticityAxiMaterial::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright,
                                                   REAL weight,
												   TPZFMatrix &ek, TPZFMatrix &ef){
	
	TPZFMatrix &dphiL = dataleft.dphix;
	TPZFMatrix &dphiR = dataright.dphix;
	TPZFMatrix &phiL = dataleft.phi;
	TPZFMatrix &phiR = dataright.phi;
	TPZFMatrix &axesL = dataleft.axes;
	TPZFMatrix &axesR = dataright.axes;
	
	TPZManVector<REAL,3> &normal = data.normal;
	
	int &LeftPOrder=dataleft.p;
	int &RightPOrder=dataright.p;
	
	REAL &faceSize=data.HSize;
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		dataleft.phi.Print("phil",sout);
		dataright.phi.Print("phir",sout);
		LOGPZ_DEBUG(logger,sout.str())
	}
#endif
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		sout << "Origin = { " << f_Origin << "};" << std::endl;
		sout << "AxisR = { " << f_AxisR << "};" << std::endl;
		sout << "AxisZ = { " << f_AxisZ << "};" << std::endl;
		data.PrintMathematica(sout);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
	int nrowl = phiL.Rows();
	int nrowr = phiR.Rows();
	int il,jl,ir,jr;	
	
	/** R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!*/
	REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];
	int s = (R > 0)? 1:-1;
	R = fabs(R);
	
	
	TPZFNMatrix<16> dphiRZiL(2,1), dphiRZjL(2,1), dphiRZiR(2,1), dphiRZjR(2,1);
	
	double axis0DOTrL = 0., axis1DOTrL = 0., axis0DOTzL = 0., axis1DOTzL = 0.;
	double axis0DOTrR = 0., axis1DOTrR = 0., axis0DOTzR = 0., axis1DOTzR = 0.;
	for(int pos = 0; pos < 3; pos++)
	{
		axis0DOTrL += axesL.GetVal(0,pos) * f_AxisR[pos] * s;
		axis1DOTrL += axesL.GetVal(1,pos) * f_AxisR[pos] * s;
		axis0DOTzL += axesL.GetVal(0,pos) * f_AxisZ[pos];
		axis1DOTzL += axesL.GetVal(1,pos) * f_AxisZ[pos];
		axis0DOTrR += axesR.GetVal(0,pos) * f_AxisR[pos] * s;
		axis1DOTrR += axesR.GetVal(1,pos) * f_AxisR[pos] * s;
		axis0DOTzR += axesR.GetVal(0,pos) * f_AxisZ[pos];
		axis1DOTzR += axesR.GetVal(1,pos) * f_AxisZ[pos];
	}
	
	REAL R2PI = 2.0 * M_PI * R;
	
	REAL symmetry = fSymmetric; //thermo of symmetry ( -1: method symmetric; 1: method not symmetric)
	REAL penalty = fPenalty; //thermo of penalty
	penalty *= (0.5 * (LeftPOrder*LeftPOrder + RightPOrder*RightPOrder)) / faceSize;
	
	REAL beta;
	if (symmetry==1.0) {
		beta=1.;
	}else {
		beta=1.0;
	}
	
	double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
	double mu =  fE/(2.*(1. + fnu));
	
    int numbersol = dataleft.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout.precision(16);
		sout << "Penalty = " << fPenalty << ";" << std::endl;
		sout << "R = " << R << ";" << std::endl;
		sout << "fE = " << fE << ";" << std::endl;
		sout << "fnu = " << fnu << ";" << std::endl;
		sout << "Lambda = " << lambda << ";" << std::endl;
		sout << "Mu = " << mu << ";" << std::endl;
		sout << "weight = " << weight << ";" << std::endl;
		data.dsol[0].Print("dsol = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
#ifdef LOG4CXX
	TPZFNMatrix<4> DSolLAxes(2,2), DSolRAxes(2,2);
	DSolLAxes(0,0) = dataleft.dsol[0](0,0)*axis0DOTrL+dataleft.dsol[0](1,0)*axis1DOTrL;
	DSolLAxes(1,0) = dataleft.dsol[0](0,0)*axis0DOTzL+dataleft.dsol[0](1,0)*axis1DOTzL;
	DSolLAxes(0,1) = dataleft.dsol[0](0,1)*axis0DOTrL+dataleft.dsol[0](1,1)*axis1DOTrL;
	DSolLAxes(1,1) = dataleft.dsol[0](0,1)*axis0DOTzL+dataleft.dsol[0](1,1)*axis1DOTzL;
	
	DSolRAxes(0,0) = dataright.dsol[0](0,0)*axis0DOTrR+dataright.dsol[0](1,0)*axis1DOTrR;
	DSolRAxes(1,0) = dataright.dsol[0](0,0)*axis0DOTzR+dataright.dsol[0](1,0)*axis1DOTzR;
	DSolRAxes(0,1) = dataright.dsol[0](0,1)*axis0DOTrR+dataright.dsol[0](1,1)*axis1DOTrR;
	DSolRAxes(1,1) = dataright.dsol[0](0,1)*axis0DOTzR+dataright.dsol[0](1,1)*axis1DOTzR;
	
	TPZFNMatrix<9> DeformL(3,3,0.),DeformR(3,3,0.);
	DeformL(0,0) = DSolLAxes(0,0); 
	DeformL(1,1) = DSolLAxes(1,1);
	DeformL(0,1) = (DSolLAxes(0,1)+DSolLAxes(1,0))/2.;
	DeformL(1,0) = DeformL(0,1);
	DeformL(2,2) = dataleft.sol[0][0]/R;
	
	DeformR(0,0) = DSolRAxes(0,0); 
	DeformR(1,1) = DSolRAxes(1,1);
	DeformR(0,1) = (DSolRAxes(0,1)+DSolRAxes(1,0))/2.;
	DeformR(1,0) = DeformR(0,1);
	DeformR(2,2) = dataright.sol[0][0]/R;
	
	TPZFNMatrix<9> TensorL(3,3,0.),TensorR(3,3,0.);
	REAL TrDeformL, TrDeformR;
	TrDeformL = DeformL(0,0)+DeformL(1,1)+DeformL(2,2);
	TrDeformR = DeformR(0,0)+DeformR(1,1)+DeformR(2,2);
	TensorL = (2.*mu)*DeformL;
	TensorR = (2.*mu)*DeformR;
	REAL normalCompL, normalCompR;
	for (int i=0; i<3; i++) {
		TensorL(i,i) += lambda*TrDeformL;
		TensorR(i,i) += lambda*TrDeformR;
	}
	normalCompL = TensorL(1,0)*normal[0]+TensorL(1,1)*normal[1];
	normalCompR = TensorR(1,0)*normal[0]+TensorR(1,1)*normal[1];
	if (logdata->isDebugEnabled() && TrDeformL != 0.)
	{
		std::stringstream sout;
		sout.precision(15);
		TensorL.Print("TensorL = ",sout,EMathematicaInput);
		TensorR.Print("TensorR = ",sout,EMathematicaInput);
		sout << "NormalCompL = " << normalCompL << ";" << std::endl;
		sout << "NormalCompR = " << normalCompR << ";" << std::endl;
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
	
	//Calcule: Integrate { -[v].<sigma(u).n> + symmetry *[u].<sigma(v).n> + penalty*[u].[v] } ds
	
	// 1) Matrix Band:  phi_I_Left, phi_J_Left
	
	for(il = 0; il< nrowl; il++ )
	{	
		//dphiL_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
		dphiRZiL.PutVal(0,0, dphiL.GetVal(0,il)*axis0DOTrL + dphiL.GetVal(1,il)*axis1DOTrL );
		
		//dphiL_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
		dphiRZiL.PutVal(1,0, dphiL.GetVal(0,il)*axis0DOTzL + dphiL.GetVal(1,il)*axis1DOTzL );
		
		
		for( jl = 0; jl < nrowl; jl++ )
		{
			//dphiL_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
			dphiRZjL.PutVal(0,0, dphiL.GetVal(0,jl)*axis0DOTrL + dphiL.GetVal(1,jl)*axis1DOTrL );
			
			//dphiL_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
			dphiRZjL.PutVal(1,0, dphiL.GetVal(0,jl)*axis0DOTzL + dphiL.GetVal(1,jl)*axis1DOTzL );
			
			double nr = normal[0];
			double nz = normal[1];
			
			double term00 = -1.*(lambda*nr*phiL(il,0)*phiL(jl,0)/(2.*R) +
								 mu*nz*phiL(il,0)*dphiRZjL(1,0)/2. + 
								 lambda*nr*phiL(il,0)*dphiRZjL(0,0)/2. + 
								 mu*nr*phiL(il,0)*dphiRZjL(0,0));
			term00 += symmetry*(lambda*nr*phiL(il,0)*phiL(jl,0)/(2.*R) +
								mu*nz*phiL(jl,0)*dphiRZiL(1,0)/2. + 
								lambda*nr*phiL(jl,0)*dphiRZiL(0,0)/2. +
								mu*nr*phiL(jl,0)*dphiRZiL(0,0));
			term00 += beta*penalty*phiL(il,0)*phiL(jl,0);
			
			
			double term01 = -1.*(lambda*nr*phiL(il,0)*dphiRZjL(1,0)/2. +
								 mu*nz*phiL(il,0)*dphiRZjL(0,0)/2.);
			term01 += symmetry*(lambda*nz*phiL(il,0)*phiL(jl,0)/(2.*R) +
								mu*nr*phiL(jl,0)*dphiRZiL(1,0)/2. + 
								lambda*nz*phiL(jl,0)*dphiRZiL(0,0)/2.);
			
			
			double term10 = -1.*(lambda*nz*phiL(il,0)*phiL(jl,0)/(2.*R) +
								 mu*nr*phiL(il,0)*dphiRZjL(1,0)/2. + 
								 lambda*nz*phiL(il,0)*dphiRZjL(0,0)/2.);
			term10 += symmetry*(lambda*nr*phiL(jl,0)*dphiRZiL(1,0)/2. + 
								mu*nz*phiL(jl,0)*dphiRZiL(0,0)/2.);
			
			
			double term11 = -1.*(lambda*nz*phiL(il,0)*dphiRZjL(1,0)/2. + 
								 mu*nz*phiL(il,0)*dphiRZjL(1,0) + 
								 mu*nr*phiL(il,0)*dphiRZjL(0,0)/2.);
			term11 += symmetry*(lambda*nz*phiL(jl,0)*dphiRZiL(1,0)/2. +
								mu*nz*phiL(jl,0)*dphiRZiL(1,0) + 
								mu*nr*phiL(jl,0)*dphiRZiL(0,0)/2.);
			term11 +=beta*penalty*phiL(il,0)*phiL(jl,0);
			
			
			ek(2*il, 2*jl)     += weight*R2PI*term00;
			ek(2*il, 2*jl+1) += weight*R2PI*term01;
			ek(2*il+1, 2*jl)     += weight*R2PI*term10;
			ek(2*il+1, 2*jl+1) += weight*R2PI*term11;
			
		}
	}
	
	
	// 2) Matrix Band:  phi_I_Left, phi_J_Right
	for( il = 0; il < nrowl; il++ )
	{	
		//dphiL_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
		dphiRZiL.PutVal(0,0, dphiL.GetVal(0,il)*axis0DOTrL + dphiL.GetVal(1,il)*axis1DOTrL );
		
		//dphiL_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
		dphiRZiL.PutVal(1,0, dphiL.GetVal(0,il)*axis0DOTzL + dphiL.GetVal(1,il)*axis1DOTzL );
		
		for(jr = 0; jr < nrowr; jr++ )
		{
			//dphiR_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
			dphiRZjR.PutVal(0,0, dphiR.GetVal(0,jr)*axis0DOTrR + dphiR.GetVal(1,jr)*axis1DOTrR );
			
			//dphiR_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
			dphiRZjR.PutVal(1,0, dphiR.GetVal(0,jr)*axis0DOTzR + dphiR.GetVal(1,jr)*axis1DOTzR );
			
			
			double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
			double mu =  fE/(2.*(1. + fnu));
			double nr = normal[0];
			double nz = normal[1];
			
			double term02 = -1.*(lambda*nr*phiL(il,0)*phiR(jr,0)/(2.*R) +
								 mu*nz*phiL(il,0)*dphiRZjR(1,0)/2. + 
								 lambda*nr*phiL(il,0)*dphiRZjR(0,0)/2. +
								 mu*nr*phiL(il,0)*dphiRZjR(0,0));
			term02 += -1.*symmetry*(lambda*nr*phiL(il,0)*phiR(jr,0)/(2.*R) +
									mu*nz*phiR(jr,0)*dphiRZiL(1,0)/2. + 
									lambda*nr*phiR(jr,0)*dphiRZiL(0,0)/2. +
									mu*nr*phiR(jr,0)*dphiRZiL(0,0));
			term02 += -1.*beta*penalty*phiL(il,0)*phiR(jr,0);
			
			
			double term03 = -1.*(lambda*nr*phiL(il,0)*dphiRZjR(1,0)/2. +
								 mu*nz*phiL(il,0)*dphiRZjR(0,0)/2.);
			term03 += -1.*symmetry*(lambda*nz*phiL(il,0)*phiR(jr,0)/(2.*R) +
									mu*nr*phiR(jr,0)*dphiRZiL(1,0)/2. + 
									lambda*nz*phiR(jr,0)*dphiRZiL(0,0)/2.);
			
			
			double term12 = -1.*(lambda*nz*phiL(il,0)*phiR(jr,0)/(2.*R) +
								 mu*nr*phiL(il,0)*dphiRZjR(1,0)/2. + 
								 lambda*nz*phiL(il,0)*dphiRZjR(0,0)/2.);
			term12 += -1.*symmetry*(lambda*nr*phiR(jr,0)*dphiRZiL(1,0)/2. +
									mu*nz*phiR(jr,0)*dphiRZiL(0,0)/2.);
			
			
			double term13 = -1.*(lambda*nz*phiL(il,0)*dphiRZjR(1,0)/2. +
								 mu*nz*phiL(il,0)*dphiRZjR(1,0) + 
								 mu*nr*phiL(il,0)*dphiRZjR(0,0)/2.);
			term13 += -1.*symmetry*(lambda*nz*phiR(jr,0)*dphiRZiL(1,0)/2. +
									mu*nz*phiR(jr,0)*dphiRZiL(1,0) + 
									mu*nr*phiR(jr,0)*dphiRZiL(0,0)/2.);
			term13 += -1.*beta*penalty*phiL(il,0)*phiR(jr,0);
			
			ek(2*il, 2*jr+2*nrowl) += weight*R2PI*term02;
			ek(2*il, 2*jr+2*nrowl+1) += weight*R2PI*term03;
			ek(2*il+1, 2*jr+2*nrowl) += weight*R2PI*term12;
			ek(2*il+1, 2*jr+2*nrowl+1) += weight*R2PI*term13;
			
		}
	}
	
	// 3) Matrix Band:  phi_I_Right, phi_J_Left
	for( ir = 0; ir < nrowr; ir++ )
	{	
		//dphiR_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
		dphiRZiR.PutVal(0,0, dphiR.GetVal(0,ir)*axis0DOTrR + dphiR.GetVal(1,ir)*axis1DOTrR );
		
		//dphiR_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
		dphiRZiR.PutVal(1,0, dphiR.GetVal(0,ir)*axis0DOTzR + dphiR.GetVal(1,ir)*axis1DOTzR );
		
		for( jl = 0; jl < nrowl; jl++ )
		{
			//dphiL_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
			dphiRZjL.PutVal(0,0, dphiL.GetVal(0,jl)*axis0DOTrL + dphiL.GetVal(1,jl)*axis1DOTrL );
			
			//dphiL_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
			dphiRZjL.PutVal(1,0, dphiL.GetVal(0,jl)*axis0DOTzL + dphiL.GetVal(1,jl)*axis1DOTzL );
			
			
			double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
			double mu =  fE/(2.*(1. + fnu));
			double nr = normal[0];
			double nz = normal[1];
			
			double term20 = (lambda*nr*phiR(ir,0)*phiL(jl,0)/(2.*R) + 
							 mu*nz*phiR(ir,0)*dphiRZjL(1,0)/2. + 
							 lambda*nr*phiR(ir,0)*dphiRZjL(0,0)/2. +
							 mu*nr*phiR(ir,0)*dphiRZjL(0,0));
			term20 += symmetry*(lambda*nr*phiR(ir,0)*phiL(jl,0)/(2.*R) +
								mu*nz*phiL(jl,0)*dphiRZiR(1,0)/2. + 
								lambda*nr*phiL(jl,0)*dphiRZiR(0,0)/2. +
								mu*nr*phiL(jl,0)*dphiRZiR(0,0));
			term20 += -1.*beta*penalty*phiR(ir,0)*phiL(jl,0);
			
			
			double term21 = (lambda*nr*phiR(ir,0)*dphiRZjL(1,0)/2. +
							 mu*nz*phiR(ir,0)*dphiRZjL(0,0)/2.);
			term21 += symmetry*(lambda*nz*phiR(ir,0)*phiL(jl,0)/(2.*R) +
								mu*nr*phiL(jl,0)*dphiRZiR(1,0)/2. + 
								lambda*nz*phiL(jl,0)*dphiRZiR(0,0)/2.);
			
			
			double term30 = (lambda*nz*phiR(ir,0)*phiL(jl,0)/(2.*R) +
							 mu*nr*phiR(ir,0)*dphiRZjL(1,0)/2. + 
							 lambda*nz*phiR(ir,0)*dphiRZjL(0,0)/2.);
			term30 += symmetry*(lambda*nr*phiL(jl,0)*dphiRZiR(1,0)/2. + 
								mu*nz*phiL(jl,0)*dphiRZiR(0,0)/2.);
			
			
			double term31 = (lambda*nz*phiR(ir,0)*dphiRZjL(1,0)/2. +
							 mu*nz*phiR(ir,0)*dphiRZjL(1,0) + 
							 mu*nr*phiR(ir,0)*dphiRZjL(0,0)/2.);
			term31 += symmetry*(lambda*nz*phiL(jl,0)*dphiRZiR(1,0)/2. +
								mu*nz*phiL(jl,0)*dphiRZiR(1,0) + 
								mu*nr*phiL(jl,0)*dphiRZiR(0,0)/2.);
			term31 += -1.*beta*penalty*phiR(ir,0)*phiL(jl,0);
			
			ek(2*ir+2*nrowl, 2*jl)     += weight*R2PI*term20;
			ek(2*ir+2*nrowl, 2*jl+1) += weight*R2PI*term21;
			ek(2*ir+2*nrowl+1, 2*jl)     += weight*R2PI*term30;
			ek(2*ir+2*nrowl+1, 2*jl+1) += weight*R2PI*term31;
			
		}
	}
	
	// 4) Matrix Band:  phi_I_Right, phi_J_Right
	
	for(ir = 0; ir < nrowr; ir++ )
	{	
		//dphiR_i/dr = dphi_i/axis0 <axes0,f_AxisR> + dphi_i/axis1 <axes1,f_AxisR>
		dphiRZiR.PutVal(0,0, dphiR.GetVal(0,ir)*axis0DOTrR + dphiR.GetVal(1,ir)*axis1DOTrR );
		
		//dphiR_i/dz = dphi_i/axis0 <axes0,f_AxisZ> + dphi_i/axis1 <axes1,f_AxisZ>
		dphiRZiR.PutVal(1,0, dphiR.GetVal(0,ir)*axis0DOTzR + dphiR.GetVal(1,ir)*axis1DOTzR );
		
		for( jr = 0; jr < nrowr; jr++ )
		{
			//dphiR_j/dr = dphi_j/axis0 <axes0,f_AxisR> + dphi_j/axis1 <axes1,f_AxisR>
			dphiRZjR.PutVal(0,0, dphiR.GetVal(0,jr)*axis0DOTrR + dphiR.GetVal(1,jr)*axis1DOTrR );
			
			//dphiR_j/dz = dphi_j/axis0 <axes0,f_AxisZ> + dphi_j/axis1 <axes1,f_AxisZ>
			dphiRZjR.PutVal(1,0, dphiR.GetVal(0,jr)*axis0DOTzR + dphiR.GetVal(1,jr)*axis1DOTzR );
			
			double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
			double mu =  fE/(2.*(1. + fnu));
			double nr = normal[0];
			double nz = normal[1];
			
			double term22 = (lambda*nr*phiR(ir,0)*phiR(jr,0)/(2.*R) +
							 mu*nz*phiR(ir,0)*dphiRZjR(1,0)/2. + 
							 lambda*nr*phiR(ir,0)*dphiRZjR(0,0)/2. +
							 mu*nr*phiR(ir,0)*dphiRZjR(0,0));
			term22 += -1.*symmetry*(lambda*nr*phiR(ir,0)*phiR(jr,0)/(2.*R) +
									mu*nz*phiR(jr,0)*dphiRZiR(1,0)/2. + 
									lambda*nr*phiR(jr,0)*dphiRZiR(0,0)/2. +
									mu*nr*phiR(jr,0)*dphiRZiR(0,0));
			term22 += beta*penalty*phiR(ir,0)*phiR(jr,0);
			
			
			double term23 = (lambda*nr*phiR(ir,0)*dphiRZjR(1,0)/2. +
							 mu*nz*phiR(ir,0)*dphiRZjR(0,0)/2.);
			term23 += -1.*symmetry*(lambda*nz*phiR(ir,0)*phiR(jr,0)/(2.*R) +
									mu*nr*phiR(jr,0)*dphiRZiR(1,0)/2. + 
									lambda*nz*phiR(jr,0)*dphiRZiR(0,0)/2.);
			
			
			double term32 = (lambda*nz*phiR(ir,0)*phiR(jr,0)/(2.*R) +
							 mu*nr*phiR(ir,0)*dphiRZjR(1,0)/2. + 
							 lambda*nz*phiR(ir,0)*dphiRZjR(0,0)/2.);
			term32 += -1.*symmetry*(lambda*nr*phiR(jr,0)*dphiRZiR(1,0)/2. +
									mu*nz*phiR(jr,0)*dphiRZiR(0,0)/2.);
			
			
			double term33 = (lambda*nz*phiR(ir,0)*dphiRZjR(1,0)/2. +
							 mu*nz*phiR(ir,0)*dphiRZjR(1,0) + 
							 mu*nr*phiR(ir,0)*dphiRZjR(0,0)/2.);
			term33 += -1.*symmetry*(lambda*nz*phiR(jr,0)*dphiRZiR(1,0)/2. +
									mu*nz*phiR(jr,0)*dphiRZiR(1,0) + 
									mu*nr*phiR(jr,0)*dphiRZiR(0,0)/2.);
			term33 += beta*penalty*phiR(ir,0)*phiR(jr,0);
			
			
			ek(2*ir+2*nrowl, 2*jr+2*nrowl) += weight*R2PI*term22;
			ek(2*ir+2*nrowl, 2*jr+2*nrowl+1) += weight*R2PI*term23;
			ek(2*ir+2*nrowl+1, 2*jr+2*nrowl) += weight*R2PI*term32;
			ek(2*ir+2*nrowl+1, 2*jr+2*nrowl+1) += weight*R2PI*term33;
			
		}
	}
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif
}

//-------------------------------------------------------------------

/** returns the variable index associated with the name*/
int TPZElasticityAxiMaterial::VariableIndex(const std::string &name)
{
	if(!strcmp("Eigenvector1",name.c_str()))     return 9;
	if(!strcmp("Eigenvector2",name.c_str()))     return 1;
	if(!strcmp("Eigenvector3",name.c_str()))     return 2;
	if(!strcmp("Sigmarr",name.c_str()))       return 3;
	if(!strcmp("Sigmazz",name.c_str()))       return 4;
	if(!strcmp("Sigmatt",name.c_str()))       return 5;
	if(!strcmp("Taurz",name.c_str()))         return 6;
	if(!strcmp("displacement",name.c_str()))  return 7;
	if(!strcmp("MohrCoulomb",name.c_str()))  return 8;
	
	return TPZMaterial::VariableIndex(name);
	return -1;
}

/**returns the number of variables associated with the
 variable indexed by var, var is obtained by calling VariableIndex*/
int TPZElasticityAxiMaterial::NSolutionVariables(int var)
{
	switch(var) {
		case 9:
			return 3;
		case 1:
			return 3;
		case 2:
			return 3;
		case 3:
			return 1;
		case 4:
			return 1;
		case 5:
			return 1;
		case 6:
			return 1;
		case 7:
			return 3;
		case 8:
			return 1;
		default:
			return TPZMaterial::NSolutionVariables(var);
			return 0;
	}
}

/** returns the solution associated with the var index based
 on the finite element approximation*/
void TPZElasticityAxiMaterial::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
	if(var == 0) 
	{
		TPZMaterial::Solution(data,var,Solout);
		return;
	}
    int numbersol = data.dsol.size();
    if (numbersol != 1) {
        DebugStop();
    }

	TPZFMatrix &axes = data.axes;
	TPZVec<REAL> &SolAxes = data.sol[0];
	TPZFMatrix &DSolAxes = data.dsol[0];
	
    // R = Dot[{data.x - origin},{AxisR}]   ***because AxisR is already normalized!
    REAL R = (data.x[0] - f_Origin[0])*f_AxisR[0] + (data.x[1] - f_Origin[1])*f_AxisR[1] + (data.x[2] - f_Origin[2])*f_AxisR[2];
    REAL Z = (data.x[0] - f_Origin[0])*f_AxisZ[0] + (data.x[1] - f_Origin[1])*f_AxisZ[1] + (data.x[2] - f_Origin[2])*f_AxisZ[2];
    
    int s = (R > 0)? 1:-1;
    R = fabs(R);
	if(R < 1.e-10) R = 1.e-10;
    
    REAL DelTemp(fDelTemperature);
    
	if (fTemperatureFunction) {
        TPZManVector<REAL,2> RZ(2);
        RZ[0] = R;
        RZ[1] = Z;
        fTemperatureFunction(RZ,DelTemp);
    }
	
	
	double axis0DOTr = 0., axis1DOTr = 0., axis0DOTz = 0., axis1DOTz = 0.;
	for(int pos = 0; pos < 3; pos++)
	{
		axis0DOTr += axes.GetVal(0,pos) * f_AxisR[pos] * s;
		axis1DOTr += axes.GetVal(1,pos) * f_AxisR[pos] * s;
		axis0DOTz += axes.GetVal(0,pos) * f_AxisZ[pos];
		axis1DOTz += axes.GetVal(1,pos) * f_AxisZ[pos];
	}
	TPZFNMatrix<9> DSolrz(2,2,0.);
	DSolrz.PutVal(0,0, DSolAxes(0,0)*axis0DOTr + DSolAxes(1,0)*axis1DOTr );
	DSolrz.PutVal(0,1, DSolAxes(0,0)*axis0DOTz + DSolAxes(1,0)*axis1DOTz );
	DSolrz.PutVal(1,0, DSolAxes(0,1)*axis0DOTr + DSolAxes(1,1)*axis1DOTr );
	DSolrz.PutVal(1,1, DSolAxes(0,1)*axis0DOTz + DSolAxes(1,1)*axis1DOTz );
	
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "Point " << data.x << std::endl;
		sout << "Solution " << data.sol << std::endl;
		DSolrz.Print("Derivatives of the solution\n",sout);
		sout << "Radius " << R << std::endl;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	//Infinitesimal Tensor
	TPZFNMatrix<9> Einf(3,3,0.);
	Einf.PutVal(0,0,DSolrz(0,0)); 
	Einf.PutVal(0,1,0.5*(DSolrz(0,1) + DSolrz(1,0)));
	Einf.PutVal(1,0,0.5*(DSolrz(0,1) + DSolrz(1,0))); 
	Einf.PutVal(1,1,DSolrz(1,1));
	Einf.PutVal(2,2,data.sol[0][0]/R);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		Einf.Print("Deformation tensor",sout);
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	double lambda = -((fE*fnu)/((1. + fnu)*(2.*fnu-1.)));
	double mi =  fE/(2.*(1. + fnu));
	double trE = Einf(0,0) + Einf(1,1) + Einf(2,2);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "E = " << fE << " fnu = " << fnu << " mi = " << mi << " lambda = " << lambda << " trE = " << trE;
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	//Stress Tensor
	TPZFNMatrix<9> T(3,3,0.);
	double cte = lambda*trE;
    REAL epsT = DelTemp*fAlpha;
    REAL TensorThermico = (3.*lambda+2.*mi)*epsT;
	
	
	//T = lambda.tr(E).I + 2.mi.E
	T.PutVal(0,0,cte + 2.*mi*Einf(0,0)-TensorThermico); 
	T.PutVal(0,1,      2.*mi*Einf(0,1)); 
	T.PutVal(0,2,      2.*mi*Einf(0,2));
	T.PutVal(1,0,      2.*mi*Einf(1,0)); 
	T.PutVal(1,1,cte + 2.*mi*Einf(1,1)-TensorThermico); 
	T.PutVal(1,2,      2.*mi*Einf(1,2));
	T.PutVal(2,0,      2.*mi*Einf(2,0)); 
	T.PutVal(2,1,      2.*mi*Einf(2,1)); 
	T.PutVal(2,2,cte + 2.*mi*Einf(2,2)-TensorThermico);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		T.Print("Stress tensor",sout);
		LOGPZ_DEBUG(logger,sout.str());
	}
#endif
	
	switch(var)
	{
		case 9: //Solout = 1stEigenvalue * {1stEigenvector}
		{
			int NumIt = 1000;
			REAL tol = 1.E-5;
			TPZVec<REAL> EigValues(3,0.);
			TPZFNMatrix<9> EigVectors(3,3,0.);
			bool EigenWorks;
			EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
			if(EigenWorks)
			{
				Solout.Resize(3);
				for(int i = 0; i < 3; i++) Solout[i] = EigValues[0] * EigVectors(0,i);
			}
			else cout << "TPZElasticityAxiMaterial::Solution Error -> case 0\n";
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "First eigenvector " << std::endl;
				sout << Solout;
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 1: //Solout = 2ndEigenvalue * {2ndEigenvector}
		{
			int NumIt = 1000;
			REAL tol = 1.E-5;
			TPZVec<REAL> EigValues(3,0.);
			TPZFNMatrix<9> EigVectors(3,3,0.);
			bool EigenWorks;
			EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
			if(EigenWorks)
			{
				Solout.Resize(3);
				for(int i = 0; i < 3; i++) Solout[i] = EigValues[1] * EigVectors(1,i);
			}
			else cout << "TPZElasticityAxiMaterial::Solution Error -> case 1\n";
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Second eigenvector " << std::endl;
				sout << Solout;
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 2: //Solout = 3rdEigenvalue * {3rdEigenvector}
		{
			int NumIt = 1000;
			REAL tol = 1.E-5;
			TPZVec<REAL> EigValues(3,0.);
			TPZFNMatrix<9> EigVectors(3,3,0.);
			bool EigenWorks;
			EigenWorks = T.SolveEigensystemJacobi(NumIt, tol, EigValues, EigVectors);
			if(EigenWorks)
			{
				Solout.Resize(3);
				for(int i = 0; i < 3; i++) Solout[i] = EigValues[2] * EigVectors(2,i);
			}
			else cout << "TPZElasticityAxiMaterial::Solution Error -> case 2\n";
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Third eigenvector " << std::endl;
				sout << Solout;
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 3: //Solout = Sigma_r,r
		{
			Solout.Resize(1);
			Solout[0] = T(0,0);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Sigma r " << T(0,0);
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 4: //Solout = Sigma_z,z
		{
			Solout.Resize(1);
			Solout[0] = T(1,1);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Sigma z " << T(1,1);
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 5: //Solout = Sigma_theta,theta
		{
			Solout.Resize(1);
			Solout[0] = T(2,2);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Sigma theta " << T(2,2);
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 6: //Solout = Tau_r,z
		{
			Solout.Resize(1);
			Solout[0] = T(0,1);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Sigma rz " << T(0,1);
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		case 7: //Solout = Displacement
		{
			Solout.Resize(3);
			Solout[0] = SolAxes[0];
			Solout[1] = SolAxes[1];
			Solout[2] = SolAxes[2];
		}
			break;
			
		case 8: //MohrCoulomb plasticity criteria
		{
			Solout.Resize(1);
			REAL i1, i2, i3, j1, j2, j3;
			
			i1 = T(0,0) + T(1,1) + T(2,2);
			TPZFMatrix T2(3,3,0.);
			T.Multiply(T,T2);
			
			i2 = 0.5*( i1*i1 - (T2(0,0) + T2(1,1) + T2(2,2)) );
			
			i3 = T(0,0)*T(1,1)*T(2,2) + T(0,1)*T(1,2)*T(2,0) + T(0,2)*T(1,0)*T(2,1) - T(0,2)*T(1,1)*T(2,0) - T(1,2)*T(2,1)*T(0,0) - T(2,2)*T(0,1)*T(1,0);
			
			j1 = 0.;
			j2 = 1./3.*(i1*i1 - 3.*i2);
			j3 = 1./27.*(2.*i1*i1*i1 - 9.*i1*i2 + 27.*i3);
			
			REAL cos3theta = 3.*sqrt(3.)/2. * j3/(pow(j2,(REAL)1.5));
			
			REAL theta = acos(cos3theta)/3.;
			
			Solout[0] = 1./3.*i1*sin(f_phi) + sqrt(i2)*sin(theta + M_PI/3.) + sqrt(i2/3.)*cos(theta + M_PI/3.)*sin(f_phi) - f_c*cos(f_phi);
#ifdef LOG4CXX
			{
				std::stringstream sout;
				sout << "Criterio de Mohr Coulomb \ni1 " << i1 << " i2 " << i2 <<
				" i3 " << i3 << " \nj1 " << j1 << " j2 " << j2 << " j3 " << j3
				<< " \nf_phi " << f_phi << 
				"  f_c " << f_c << " theta " << theta << " MohrCoul " << Solout[0];
				LOGPZ_DEBUG(logger,sout.str());
			}
#endif
		}
			break;
			
		default:
		{
			cout << "TPZElasticityAxiMaterial::Solution Error -> default\n";
			TPZMaterial::Solution(data.sol[0],data.dsol[0],data.axes,var,Solout);
		}
			break;
	}
}

void TPZElasticityAxiMaterial::Flux(TPZVec<REAL> &x, TPZVec<REAL> &Sol, TPZFMatrix &DSol, TPZFMatrix &axes, TPZVec<REAL> &flux)
{
	if(fabs(axes(2,0)) >= 1.e-6 || fabs(axes(2,1)) >= 1.e-6)
	{
		cout << "TPZElasticityAxiMaterial::Flux only serves for xy configuration\n";
		axes.Print("axes");
	}
}

void TPZElasticityAxiMaterial::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudaxes, TPZFMatrix &axes, TPZVec<REAL> &flux,
									  TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values)
{
	values[0] = 0.;
	TPZManVector<REAL> sigma(3,0.),sigma_exact(3,0.);
	REAL sigx,sigy,sigxy,gamma;
	
	TPZFNMatrix<9> du;//du = dudx
	TPZAxesTools::Axes2XYZ(dudaxes, du, axes);
	
	//tensoes aproximadas : uma forma
	gamma = du(1,0)+du(0,1);
	sigma[0] = fEover1MinNu2*(du(0,0)+fnu*du(1,1));
	sigma[1] = fEover1MinNu2*(fnu*du(0,0)+du(1,1));
	sigma[2] = fE*0.5/(1.+fnu)*gamma;

	TPZMaterialData mydata;
	mydata.sol[0]  = u;
	mydata.dsol[0] = du;
	mydata.axes = axes;
	mydata.x = x;
	
	//tensoes aproximadas : outra forma
	TPZVec<REAL> sol(1);
	Solution(mydata,5,sol);
	sigma[0] = sol[0];
	Solution(mydata,6,sol);
	sigma[1] = sol[0];
	Solution(mydata,8,sol);
	sigma[2] = sol[0];
	
	//exata
	gamma = du_exact(1,0)+du_exact(0,1);
	sigma_exact[0] = fEover1MinNu2*(du_exact(0,0)+fnu*du_exact(1,1));
	sigma_exact[1] = fEover1MinNu2*(fnu*du_exact(0,0)+du_exact(1,1));
	sigma_exact[2] = fE*0.5/(1.+fnu)*gamma;
	sigx  = (sigma[0] - sigma_exact[0]);
	sigy  = (sigma[1] - sigma_exact[1]);
	sigxy = (sigma[2] - sigma_exact[2]);
	//values[0] = calculo do erro estimado em norma Energia
	values[0] = fE*(sigx*sigx + sigy*sigy + 2*fnu*sigx*sigy)/(1-fnu*fnu);
	values[0] = (values[0] + .5*fE*sigxy*sigxy/(1+fnu));
	
	//values[1] : erro em norma L2 em tensoes
	//values[1] = sigx*sigx + sigy*sigy + sigxy*sigxy;
	
	//values[1] : erro em norma L2 em deslocamentos
	values[1] = pow(fabs(u[0] - u_exact[0]),(REAL)2.0)+pow(fabs(u[1] - u_exact[1]),(REAL)2.0);  // It is important when REAL is long double, because fabs(...) returns long double then pow() must to return long double - Jorge
	
	//values[2] : erro estimado
	values[2] = 0.;
}


int TPZElasticityAxiMaterial::ClassId() const
{
	return TPZELASTICITYMATERIALID;
}

#ifndef BORLAND
template class TPZRestoreClass<TPZElasticityAxiMaterial,TPZELASTICITYMATERIALID>;
#endif

void TPZElasticityAxiMaterial::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf,context);
	buf.Read(&fE,1);
	buf.Read(&fnu,1);
	buf.Read(&fEover21PlusNu,1);
	buf.Read(&fEover1MinNu2,1);
	
	buf.Read(ff,3);
}

void TPZElasticityAxiMaterial::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf,withclassid);
	buf.Write(&fE,1);
	buf.Write(&fnu,1);
	buf.Write(&fEover21PlusNu,1);
	buf.Write(&fEover1MinNu2,1);
	
	buf.Write(ff,3);
}
