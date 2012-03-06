//$Id: pzporous.cpp,v 1.13 2010-06-11 22:13:02 diogo Exp $

#include "pzporous.h"
#include "pzmaterialid.h"
#include "poroelastoplasticid.h"
#include "pzbndcond.h"
#include "TPZLadeKim.h"  
#include "TPZSandlerDimaggio.h"
#include "pzelastoplastic.h"
#include "TPZYCDruckerPrager.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"


#ifdef LOG4CXX
#include "pzlog.h"
static LoggerPtr porousLogger(Logger::getLogger("material.pzPoro"));
#endif


template <class T, class TMEM>
TPZMatPorous<T, TMEM >::TPZMatPorous() : TBASEPOROUS(T, TMEM)(), fk(0.), fMu(1.), fStorageEps(0.), fAlpha(1.), fRhof(0.)
{
	fDeltaT = 1.;
	fTime = Advanced_CT;	
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZMatPorous<TBASEPOROUS(T, TMEM)>() constructor called ***";
    LOGPZ_INFO(porousLogger,sout.str().c_str());
  }
#endif
	
}

template <class T, class TMEM>
TPZMatPorous<T, TMEM >::TPZMatPorous(int id) : TBASEPOROUS(T, TMEM)(id), fk(0.), fMu(1.), fStorageEps(0.), fAlpha(1.), fRhof(0.)
{
	fDeltaT = 1.;
	fTime = Advanced_CT;
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZMatPorous<TBASEPOROUS(T, TMEM)>(int id) constructor called with id = " << id << " ***";
    LOGPZ_INFO(porousLogger,sout.str().c_str());
  }
#endif
	
}

template <class T, class TMEM>
TPZMatPorous<T, TMEM >::TPZMatPorous(const TPZMatPorous<T, TMEM> &mat) : TBASEPOROUS(T, TMEM)(mat), 
                               fk(mat.fk), fMu(mat.fMu),
                               fStorageEps(mat.fStorageEps), fAlpha(mat.fAlpha),
 							   fRhof(mat.fRhof)
{
	fDeltaT = mat.fDeltaT;
	fTime = mat.fTime;
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZMatPorous<T>() copy constructor called ***";
    LOGPZ_INFO(porousLogger,sout.str().c_str());
  }
#endif
}

template <class T, class TMEM>
TPZMatPorous<T, TMEM >::~TPZMatPorous()
{

}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Print(std::ostream &out, const int memory)
{
	out << this->Name();
	out << "\n with template argurment TBASEPOROUS(T, TMEM) = " << TBASEPOROUS(T, TMEM)::Name();
	out << "\n Delta Time: " << fDeltaT;
	out << "\n Permeability: " << fk;
	out << "\n Fluid viscosity: " << fMu;
	out << "\n Porous medium constant strain Storage Coeff: " << fStorageEps;
	out << "\n Biot-Willis alpha: " << fAlpha;
	out << "\n Base material Data:\n";
	TBASEPOROUS(T, TMEM)::Print(out, memory);
}

template <class T, class TMEM>
int TPZMatPorous<T, TMEM >::VariableIndex(const std::string &name)
{
   if(!strcmp("PorePressure", name.c_str()))return TPZMatPorous<T, TMEM >::EPorePressure;
   return TBASEPOROUS(T, TMEM)::VariableIndex(name);
}

template <class T, class TMEM>
int TPZMatPorous<T, TMEM >::NSolutionVariables(int var)
{
   if(var == TPZMatPorous<T, TMEM >::EPorePressure)return 1;

   return TBASEPOROUS(T, TMEM)::NSolutionVariables(var);
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Solution(TPZMaterialData &data, int var, TPZVec<REAL> &Solout)
{
  if(var == TPZMatPorous<T, TMEM >::EPorePressure)
	{
		REAL Pp;
		TPZVec<REAL> dPp;
		ComputePorePressure(data, Pp, dPp);
		Solout[0] = Pp;
		return;
	}
  
  TBASEPOROUS(T, TMEM)::Solution(data, var, Solout);
	
/*
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "<<< TPZMatPorous<T>::Solution() *** Sol = " << Solout;
	
    LOGPZ_DEBUG(porousLogger,sout.str().c_str());
  }
#endif
*/
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ek, TPZFMatrix &ef)
{
  int in, jn;
  REAL val, Pp;
  TPZVec<REAL> dPp;
	
  TPZFMatrix &dphi = data.dphix, dphiXYZ;
  TPZFMatrix &phi  = data.phi;
  TPZFMatrix &axes = data.axes, axesT;
  TPZManVector<REAL,3> &x = data.x;

  const int phr = phi.Rows();
  if(this->fForcingFunction)
     this->fForcingFunction->Execute(x,this->fForce);
    
  int dim = Dimension();
  int nstate = NStateVariables();
	
  ComputePorePressure(data, Pp, dPp);
	
  if(TBASEPOROUS(T, TMEM)::fUpdateMem)
	  UpdatePorePressure(data);

	
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZMatPorous<T, TMEM >::Contribute ***";
	if(fTime == Last_CT) sout << " Last State Contribution";
	if(fTime == Advanced_CT) sout << " Advanced State Contribution";
 	LOGPZ_DEBUG(porousLogger,sout.str().c_str());
  }
#endif

//	cout << "\nPp = " << Pp;
	
  // rotating the shape functions to the XYZ coordinates
  axes.Transpose(&axesT);
  axesT.Multiply(dphi,dphiXYZ);	
  	
  TPZManVector<REAL, 3> Q(3,0);

  for(in = 0; in < phr; in++) { //in: test function index
	
	//qh contribution
	// fForce represents the gravity acceleration
	val = fRhof * (fk / fMu) * 
		  ( TBASEPOROUS(T, TMEM)::fForce[0] * dphiXYZ(0,in) +
		    TBASEPOROUS(T, TMEM)::fForce[1] * dphiXYZ(1,in) +
			TBASEPOROUS(T, TMEM)::fForce[2] * dphiXYZ(2,in) )
		  * fDeltaT;
	//qS contribution (referring to the deltaP iterative solution)
	val -= fStorageEps * phi(in,0) * data.sol[0][dim];// / fDeltaT;
	//qH contribution 
	val -= (fk / fMu) * 
		   ( dphiXYZ(0,in)*dPp[0] +
		     dphiXYZ(1,in)*dPp[1] +
		     dphiXYZ(2,in)*dPp[2] )
		  * fDeltaT;
	//qQT  (referring to the deltaP iterative solution)
	val -= fAlpha *
		   ( phi(in, 0) * data.dsol[0](0, 0) +
	         phi(in, 0) * data.dsol[0](1, 1) +
	         phi(in, 0) * data.dsol[0](2, 2) );// /fDeltaT
	  
	ef(in*nstate+dim,0) += weight * val;
	  
	//fq contributions
	Q[0] = fAlpha * Pp * dphiXYZ(0, in);
	Q[1] = fAlpha * Pp * dphiXYZ(1, in);
	Q[2] = fAlpha * Pp * dphiXYZ(2, in);
	  
	ef(in*nstate+0,0) += weight * Q[0];
	ef(in*nstate+1,0) += weight * Q[1];
	ef(in*nstate+2,0) += weight * Q[2]; 
	  
    for( jn = 0; jn < phr; jn++ ) { //jn: trial function index
 
	  // -Q matrix contributions		
	  Q[0] = fAlpha * weight * phi(jn, 0) * dphiXYZ(0, in);
	  Q[1] = fAlpha * weight * phi(jn, 0) * dphiXYZ(1, in);
	  Q[2] = fAlpha * weight * phi(jn, 0) * dphiXYZ(2, in);
	
	  ek(in * nstate + 0, jn * nstate + dim) -= Q[0];  
	  ek(in * nstate + 1, jn * nstate + dim) -= Q[1];
	  ek(in * nstate + 2, jn * nstate + dim) -= Q[2];

	  // Transpose[Q] matrix contributions
	  Q[0] = fAlpha * weight * phi(in, 0) * dphiXYZ(0, jn);
	  Q[1] = fAlpha * weight * phi(in, 0) * dphiXYZ(1, jn);
	  Q[2] = fAlpha * weight * phi(in, 0) * dphiXYZ(2, jn);	
		
	  ek(in * nstate + dim, jn * nstate + 0) += Q[0];// / fDeltaT;  
	  ek(in * nstate + dim, jn * nstate + 1) += Q[1];// / fDeltaT;
	  ek(in * nstate + dim, jn * nstate + 2) += Q[2];// / fDeltaT;

	  // S matrix contibution
	  val = fStorageEps * phi(in,0) * phi(jn,0);// / fDeltaT;
		
	  // H matrix contributions
	  val += (fk / fMu) *
			 ( dphiXYZ(0,in)*dphiXYZ(0,jn) +
			   dphiXYZ(1,in)*dphiXYZ(1,jn) +
			   dphiXYZ(2,in)*dphiXYZ(2,jn) )
			 * fDeltaT;

	  ek(in * nstate + dim, jn * nstate + dim) += weight * val;
      
    }//jn
	  
  }//in
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "*** TPZMatPorous<T>::Contribute ***";
		sout << "ek Matrix before base classe contribution: ";
		//ek.Print("",sout,EMatlabNonZeros);
		LOGPZ_DEBUG(porousLogger,sout.str().c_str());
	}
#endif
	
	TBASEPOROUS(T, TMEM)::Contribute(data, weight, ek, ef);
	
#ifdef LOG4CXX
	{
		std::stringstream sout;
		sout << "<<< TPZMatPorous<T>::Contribute ***";
		sout << "ek Matrix after base classe contribution: ";
		ek.Print("",sout,EMatlabNonZeros);
		LOGPZ_DEBUG(porousLogger,sout.str().c_str());
	}
#endif
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::ContributeBC(TPZMaterialData &data,
				                       REAL weight,
									   TPZFMatrix &ek,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
{
	
#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << ">>> TPZMatPorous<T>::ContributeBC *** with bc.Type()=" << bc.Type();
    LOGPZ_DEBUG(porousLogger,sout.str().c_str());
  }
#endif
  TPZFMatrix &phi = data.phi;

  const REAL BIGNUMBER  = 1.e12;
	
  int dim = Dimension();
  int nstate = NStateVariables();

  const int phr = phi.Rows();
  int in,jn;
  REAL v2, v1;
	
  v1 = bc.Val1()(dim, dim);
  v2 = bc.Val2()(dim, 0);
	
  switch (bc.Type()) {
  case 10: // Dirichlet Pressure condition
  case 3: // Mechanical Directional Dirichlet condition
    for(in = 0 ; in < phr; in++) {
      ef(nstate * in + dim,0) += BIGNUMBER * (v1 - v1) * phi(in,0) * weight * v2;     
      for (jn = 0 ; jn < phr; jn++) {
        ek(nstate * in + dim,nstate * jn + dim) += BIGNUMBER * phi(in,0) * phi(jn,0) * weight * v2;
      }//jn
    }//in
    break;

  case 11: // Neumann condition (boundary flow condition) (qGamma)
    for(in = 0 ; in < phi.Rows(); in++) 
      ef(nstate * in + dim,0) += - v1 * phi(in,0) * weight;
    break;
  default:
		  break;

#ifdef LOG4CXX
  {
    std::stringstream sout;
    sout << "<<< TPZMatPorous<T>::ContributeBC *** No Flow BC of Type " << bc.Type()
		 << " - Verifying mechanical BC Types in the parent class...";
    LOGPZ_DEBUG(porousLogger,sout.str().c_str());
  }
#endif

  }//switch

  TBASEPOROUS(T, TMEM)::ContributeBC(data, weight, ek, ef, bc);

}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix &ef)
{
	TPZMaterial::Contribute(data, weight, ef);//not efficient but here to remember reimplementing it when Contribute becomes robust
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::ContributeBC(TPZMaterialData &data,
									   REAL weight,
									   TPZFMatrix &ef,
									   TPZBndCond &bc)
{
    TPZMaterial::ContributeBC(data, weight, ef, bc);//not efficient but here to remember reimplementing it when ContributeBC becomes robust 
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Errors(TPZVec<REAL> &x,TPZVec<REAL> &u, TPZFMatrix &dudx, 
                    TPZFMatrix &axes, TPZVec<REAL> &flux,
                    TPZVec<REAL> &u_exact,TPZFMatrix &du_exact,TPZVec<REAL> &values)
{
	TBASEPOROUS(T, TMEM)::Errors(x,u,dudx, axes, flux, u_exact, du_exact, values);
                          
}

	
template <class T, class TMEM>
TPZAutoPointer<TPZMaterial> TPZMatPorous<T, TMEM >::NewMaterial()
{
	return new TPZMatPorous<T, TMEM>(*this);
}

template <class T, class TMEM>
int TPZMatPorous<T, TMEM >::ClassId() const
{
	return TBASEPOROUS(T, TMEM)::ClassId() - NUMPLASTICMODELS;
	// allowing different IDs for each template instantiation.
}

template <class T, class TMEM>
std::string TPZMatPorous<T, TMEM >::Name()
{
	return "TPZMatPorous<TBASEPOROUS(T, TMEM)>"; 
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Write(TPZStream &buf, int withclassid)
{
	//this->TPZSaveable::Write(buf, withclassid);

    TBASEPOROUS(T, TMEM)::Write(buf, 0);

    buf. Write(&fDeltaT, 1);
	buf. Write(&fk, 1);
	buf. Write(&fMu, 1);
	buf. Write(&fStorageEps, 1);
	buf. Write(&fAlpha, 1);

}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::Read(TPZStream &buf, void *context)
{
   // this->TPZSaveable::Read(buf, context);
	
	TBASEPOROUS(T, TMEM)::Read(buf, context);
	
    buf. Read(&fDeltaT, 1);
	buf. Read(&fk, 1);
	buf. Read(&fMu, 1);
	buf. Read(&fStorageEps, 1);
	buf. Read(&fAlpha, 1);
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::SetUp(const REAL &k, const REAL &Mu, 
								const REAL &StorageEps,
								const REAL &Alpha,
								const REAL &Rhof)
{
	fk      = k;
	fMu     = Mu;
	fStorageEps = StorageEps;
	fAlpha  = Alpha;
	fRhof   = Rhof;
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::FillDataRequirements(TPZMaterialData &data){
  	
	TBASEPOROUS(T, TMEM)::FillDataRequirements(data);
	
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::ComputePorePressure(TPZMaterialData & data, REAL & Pp, TPZVec<REAL> & dPp)
{
	int dim = Dimension(), i;
	int intPt = data.intPtIndex;
	
	// Retrieving information at time n
	Pp = TBASEPOROUS(T, TMEM)::fMemory[intPt].fPorePressure;
	dPp = TBASEPOROUS(T, TMEM)::fMemory[intPt].fdPorePressure;
	
	// adding deltaP information from n+1 time
	Pp += data.sol[0][dim];
	for(i = 0; i < dim; i++)dPp[i] += data.dsol[0](i, dim);
}

template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::UpdatePorePressure(TPZMaterialData & data)
{
	int dim = Dimension(), i;
	int intPt = data.intPtIndex;
	
	// updating n+1 information
	TBASEPOROUS(T, TMEM)::fMemory[intPt].fPorePressure += data.sol[0][dim];
		
    for(i = 0; i < dim; i++)
	{
		TBASEPOROUS(T, TMEM)::fMemory[intPt].fdPorePressure[i] += data.dsol[0](i, dim); 
	}
}


template <class T, class TMEM>
void TPZMatPorous<T, TMEM >::SetPorePressure(const REAL Pp)
{
	TBASEPOROUS(T, TMEM)::fDefaultMem.fPorePressure = Pp;
}

/*
template class TPZMatPorous< TPZMatElastoPlastic<TPZLadeKim, TPZPoroElastoPlasticMem> >;
template class TPZMatPorous< TPZMatElastoPlastic<TPZSandlerDimaggio, TPZPoroElastoPlasticMem> >;
*/

template class TPZMatPorous< TPZLadeKim, TPZPoroElastoPlasticMem >;
template class TPZMatPorous< TPZSandlerDimaggio, TPZPoroElastoPlasticMem >;
template class TPZMatPorous<TPZPlasticStep<TPZYCDruckerPrager, TPZThermoForceA, TPZElasticResponse> , TPZPoroElastoPlasticMem>;
