//
//  pzmultiphase.h
//  PZ
//
//  Created by Omar Duran on 19/08/2013.
//  Copyright (c) 2012 LabMec-Unicamp. All rights reserved.
//

#include "pzlog.h"
#include "pzmultiphase.h"
#include "pzbndcond.h"
#include "pzfmatrix.h"
#include "pzaxestools.h"


#include <iostream>

#ifdef LOG4CXX
static LoggerPtr logger(Logger::getLogger("pz.multiphase"));
#endif

#ifdef LOG4CXX
static LoggerPtr logdata(Logger::getLogger("pz.material.multiphase.data"));
#endif

TPZMultiphase::TPZMultiphase(): TPZDiscontinuousGalerkin()
{
    fDim = 2;
	fTheta = 0.5;	
    fDeltaT = 0.;
    fmatId = 0;
}

TPZMultiphase::TPZMultiphase(int matid, int dim): TPZDiscontinuousGalerkin(matid)
{
	// Two-dimensional analysis
	
    if(dim<0 || dim >2){
        DebugStop();
    }
    
    fDim = dim;
	fTheta = 0.5;
    fDeltaT = 0.;
    fmatId = matid;
}


TPZMultiphase::~TPZMultiphase()
{
}

int TPZMultiphase::Dimension() {return fDim;}

int TPZMultiphase::MatId() {return fmatId;}

int TPZMultiphase::NStateVariables() {return 3;}

void TPZMultiphase::Print(std::ostream &out) {
	out << "name of material : " << Name() << "\n";
	out << "Coeficient which multiplies the gradient operator "<< "my var" << std::endl;
	out << "Base Class properties :";
	TPZMaterial::Print(out);
	out << "\n";
}



// Data set

/** Capilar pressure \f$ pc = pc( Sw ) \f$ */
void TPZMultiphase::CapillaryPressure(REAL So, REAL &pc, REAL &DpcDSo){
	pc = 0.0;
	DpcDSo = 0.0;
}

/** Oil relative permeability \f$ Kro = 1 - Sw \f$ */
void TPZMultiphase::Kro(REAL Sw, REAL &Kro, REAL &dKroDSw){
	
	// linear model  
	Kro = 1.0-Sw;
	dKroDSw = -1.0;
	
	// Non-linear model
	Kro = (1-Sw)*(1-Sw);
	dKroDSw = 2.0*(1-Sw)*(-1.0);
}


/** Water relative permeability \f$ Krw = Sw \f$ */
void TPZMultiphase::Krw(REAL Sw, REAL &Krw, REAL &dKrwDSw){
  
	// linear model    
	Krw = Sw;
	dKrwDSw = 1.0;
	
	// Non-linear model
	Krw = (Sw)*(Sw);
	dKrwDSw = 2.0*(Sw)*(1.0);	
	
}


/** Rock porosity \f$ Phi = Phi( P )  */
void TPZMultiphase::Porosity(REAL po, REAL &poros, REAL &dPorosDpo){
	const REAL comp = 0.0e-10;
	const REAL pref = 1.0e6;
	const REAL porosRef = 0.1;
	poros = porosRef*exp(comp*(po-pref));
	dPorosDpo = comp*porosRef*exp(comp*(po-pref));
}


/** Oil density  \f$ RhoOil = RhoOil( P )  */
void TPZMultiphase::RhoOil(REAL po, REAL &RhoOil, REAL &dRhoOilDpo){
	const REAL Oilcomp = 0.0e-8;
	const REAL pref = 1.0e6;
	RhoOil = RhoOilSC()*exp(Oilcomp*(po-pref));
	dRhoOilDpo = Oilcomp*RhoOilSC()*exp(Oilcomp*(po-pref));
}

/** Water density  \f$ RhoWater = RhoWater( P )  */
void TPZMultiphase::RhoWater(REAL po, REAL &RhoWater, REAL &dRhoWaterDpo){
	const REAL Watercomp = 0.0e-9;
	const REAL pref = 1.0e6;
	RhoWater = RhoWaterSC()*exp(Watercomp*(po-pref));
	dRhoWaterDpo = Watercomp*RhoWaterSC()*exp(Watercomp*(po-pref));
}


/** Oil viscosity  \f$ OilViscosity = OilViscosity( P )  */
void TPZMultiphase::OilViscosity(REAL po, REAL &OilViscosity, REAL &dOilViscosityDpo){
	const REAL OilViscRef = 1.0e-3;
	OilViscosity = OilViscRef;
	dOilViscosityDpo = 0;
}

/** Water viscosity  \f$ WaterViscosity = WaterViscosity( P )  */
void TPZMultiphase::WaterViscosity(REAL po, REAL &WaterViscosity, REAL &dWaterViscosityDpo){
	const REAL WaterViscRef = 1.0e-3;
	WaterViscosity = WaterViscRef;
	dWaterViscosityDpo = 0;
}

/** Oil mobility. \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$  */
void TPZMultiphase::OilLabmda(REAL &OilLabmda, REAL Po, REAL Sw, REAL &dOilLabmdaDPo, REAL &dOilLabmdaDSw){
	
	REAL krOil,Oilviscosity,OilDensity;
	REAL dKroDSw,DOilviscosityDp,DOilDensityDp;
	
	Kro(Sw, krOil, dKroDSw);
	OilViscosity(Po, Oilviscosity,DOilviscosityDp);
	RhoOil(Po, OilDensity,DOilDensityDp);
	
	OilLabmda = ((krOil)/(Oilviscosity))*(OilDensity);
	dOilLabmdaDPo = (DOilDensityDp/Oilviscosity)*(krOil);
	dOilLabmdaDSw = (OilDensity/Oilviscosity)*(dKroDSw);
}


/** Water mobility. \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::WaterLabmda(REAL &WaterLabmda, REAL Pw, REAL Sw, REAL &dWaterLabmdaDPw, REAL &dWaterLabmdaDSw){
	
	REAL krWater,Waterviscosity,WaterDensity;
	REAL dKrwDSw,DWaterviscosityDp,DWaterDensityDp;
	
	Krw(Sw, krWater, dKrwDSw);
	WaterViscosity(Pw, Waterviscosity,DWaterviscosityDp);
	RhoWater(Pw, WaterDensity,DWaterDensityDp);
	
	WaterLabmda = ((krWater)/(Waterviscosity))*(WaterDensity);
	dWaterLabmdaDPw = (DWaterDensityDp/Waterviscosity)*(krWater);
	dWaterLabmdaDSw = (WaterDensity/Waterviscosity)*(dKrwDSw);
}



/** Bulk mobility. \f$ \lambda = \lambda( pw , Sw ) \f$  */
void TPZMultiphase::Labmda(REAL &Labmda, REAL Pw, REAL Sw, REAL &dLabmdaDPw, REAL &dLabmdaDSw){
	
	REAL OilLabmdav;
	REAL dOilLabmdaDpo,dOilLabmdaDSw;
	
	REAL WaterLabmdav;
	REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
	
	OilLabmda(OilLabmdav, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
	WaterLabmda(WaterLabmdav, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);	
		
	Labmda = OilLabmdav + WaterLabmdav;
	dLabmdaDPw = dOilLabmdaDpo+dWaterLabmdaDpo;
	dLabmdaDSw = dOilLabmdaDSw+dWaterLabmdaDSw;
}

/** Oil fractional flux. \f$ f_{Oil} = f_{Oil}( pw , Sw ) \f$  */
void TPZMultiphase::fOil(REAL &fOil, REAL Pw, REAL Sw, REAL &dfOilDPw, REAL &dfOilDSw){
	
	REAL OilLabmdab;
	REAL dOilLabmdaDpo,dOilLabmdaDSw;
	
// 	REAL WaterLabmdav;
// 	REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
	
	REAL Lambdab;
	REAL dLabmdaDPw, dLabmdaDSw;
	
	OilLabmda(OilLabmdab, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
// 	WaterLabmda(WaterLabmdav, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);
	Labmda(Lambdab,Pw,Sw,dLabmdaDPw,dLabmdaDSw);
	
	fOil = ((OilLabmdab)/(Lambdab));
	dfOilDPw = ((dOilLabmdaDpo)/(Lambdab))-((OilLabmdab)/(Lambdab*Lambdab))*dLabmdaDPw;
	dfOilDSw = ((dOilLabmdaDSw)/(Lambdab))-((OilLabmdab)/(Lambdab*Lambdab))*dLabmdaDSw;
	
}


/** Water fractional flux. \f$ f_{Water} = f_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::fWater(REAL &fWater, REAL Pw, REAL Sw, REAL &dfWaterDPw, REAL &dfWaterDSw){
	
// 	REAL OilLabmdab;
// 	REAL dOilLabmdaDpo,dOilLabmdaDSw;
	
	REAL WaterLabmdab;
	REAL dWaterLabmdaDpo,dWaterLabmdaDSw;
	
	REAL Lambdab;
	REAL dLabmdaDPw, dLabmdaDSw;
	
// 	OilLabmda(OilLabmdab, Pw, Sw, dOilLabmdaDpo, dOilLabmdaDSw);
	WaterLabmda(WaterLabmdab, Pw, Sw, dWaterLabmdaDpo, dWaterLabmdaDSw);
	Labmda(Lambdab,Pw,Sw,dLabmdaDPw,dLabmdaDSw);
		   
	fWater = ((WaterLabmdab)/(Lambdab));
	dfWaterDPw = ((dWaterLabmdaDpo)/(Lambdab))-((WaterLabmdab)/(Lambdab*Lambdab))*dLabmdaDPw;
	dfWaterDSw = ((dWaterLabmdaDSw)/(Lambdab))-((WaterLabmdab)/(Lambdab*Lambdab))*dLabmdaDSw;

}


// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////

#ifdef _AUTODIFF

/** Capilar pressure \f$ pc = pc( Sw ) \f$ */
void TPZMultiphase::CapillaryPressure(BFadREAL So, BFadREAL &pc){
	pc = 0.0;
}

/** Oil relative permeability \f$ Kro = 1 - Sw \f$ */
void TPZMultiphase::Kro(BFadREAL Sw, BFadREAL &Kro){
	Kro = 1.0-Sw;
}


/** Water relative permeability \f$ Krw = Sw \f$ */
void TPZMultiphase::Krw(BFadREAL Sw, BFadREAL &Krw){
	Krw = Sw;
}


/** Rock porosity \f$ Phi = Phi( P )  */
void TPZMultiphase::Porosity(BFadREAL po, BFadREAL &poros){
	const REAL comp = 1.0e-10;
	const REAL pref = 1.0e6;
	const REAL porosRef = 0.3;
	poros = porosRef*exp(comp*((po.val())-pref));
}


/** Oil density  \f$ RhoOil = RhoOil( P )  */
void TPZMultiphase::RhoOil(BFadREAL po, BFadREAL &RhoOil){
	const REAL Oilcomp = 0.0e-8;
	const REAL pref = 1.0e6;
	RhoOil = RhoOilSC()*exp(Oilcomp*((po.val())-pref));
}

/** Water density  \f$ RhoWater = RhoWater( P )  */
void TPZMultiphase::RhoWater(BFadREAL po, BFadREAL &RhoWater){
	const REAL Watercomp = 0.0e-9;
	const REAL pref = 1.0e6;
	RhoWater = RhoWaterSC()*exp(Watercomp*((po.val())-pref));
}


/** Oil viscosity  \f$ OilViscosity = OilViscosity( P )  */
void TPZMultiphase::OilViscosity(BFadREAL po, BFadREAL &OilViscosity){
	const REAL OilViscRef = 1.0e-3;
	OilViscosity = OilViscRef;
}

/** Water viscosity  \f$ WaterViscosity = WaterViscosity( P )  */
void TPZMultiphase::WaterViscosity(BFadREAL po, BFadREAL &WaterViscosity){
	const REAL WaterViscRef = 1.0e-3;
	WaterViscosity = WaterViscRef;
}

/** Oil mobility. \f$ \lambda_{Oil} = \lambda_{Oil}( po , Sw ) \f$  */
void TPZMultiphase::OilLabmda(BFadREAL OilLabmda, BFadREAL Po, BFadREAL &Sw){
	
	BFadREAL krOil,Oilviscosity,OilDensity;
	
	Kro(Sw, krOil);
	OilViscosity(Po,Oilviscosity);
	RhoOil(Po, OilDensity);
	
	OilLabmda = ((krOil)/(Oilviscosity))*(OilDensity);
	
}


/** Water mobility. \lambda_{Water} = \lambda_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::WaterLabmda(BFadREAL WaterLabmda, BFadREAL Pw, BFadREAL &Sw){
	
	BFadREAL krWater,Waterviscosity,WaterDensity;
	
	Krw(Sw, krWater);
	WaterViscosity(Pw,Waterviscosity);
	RhoWater(Pw,WaterDensity);
	
	WaterLabmda = ((krWater)/(Waterviscosity))*(WaterDensity);
	
}


/** Bulk mobility. \f$ \lambda = \lambda( pw , Sw ) \f$  */
void TPZMultiphase::Labmda(BFadREAL Labmda, BFadREAL Pw, BFadREAL &Sw){
	
	BFadREAL OilLabmdaf, WaterLabmdaf;
	
	OilLabmda(OilLabmdaf, Pw, Sw);
	WaterLabmda(WaterLabmdaf, Pw, Sw);	
	
	Labmda = OilLabmdaf + WaterLabmdaf;
	
}


/** Oil fractional flux. \f$ f_{Oil} = f_{Oil}( pw , Sw ) \f$  */
void TPZMultiphase::fOil(BFadREAL fOil, BFadREAL Pw, BFadREAL &Sw)
{
}


/** Water fractional flux. \f$ f_{Water} = f_{Water}( pw , Sw ) \f$  */
void TPZMultiphase::fWater(BFadREAL fWater, BFadREAL Pw, BFadREAL &Sw)
{
		
}

#endif 

// Fad Methods ///////////////////////////////////////////////////////////////////////////////////////////////////////


// Constant data

/** Oil density at standar conditions - kg/m3 */
REAL TPZMultiphase::RhoOilSC(){
	return 1000.0;
}

/** Water density at standar conditions - kg/m3 */
REAL TPZMultiphase::RhoWaterSC(){
	return 1000.0;
}

/** Gravity */
REAL TPZMultiphase::g(){
	return 0.0;
}

/** Permeabilidade absoluta */
void TPZMultiphase::K(TPZFMatrix<REAL> &K){
	K.Resize(3,3);
	K.Zero();
	K(0,0) = 1.0e-13;
	K(1,1) = 1.0e-13;
	K(2,2) = 0.0e-13;
	
}


/** @brief Absolute permeability. */
void TPZMultiphase::LoadKMap(std::string MaptoRead)
{

	std::string FileName;
	std::string stringTemp;
	FileName = MaptoRead;
	
	// Definitions
	long numelements=0;
	long numKData=0;	
	
	//	Scanning for total Number geometric elements
	long NumEntitiestoRead;
	TPZStack <std::string> SentinelString;
	{
		
		// reading a general mesh information by filter
		std::ifstream read (FileName.c_str());
		std::string FlagString;
		int flag = -1;		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			if(str == "KABSOLUTE" || str == "KABSOLUTE\r") 
			{
				SentinelString.push_back(str);
			}			
			
		}
		
		FlagString = "EndReading";
		SentinelString.push_back(FlagString);
	}
	
	NumEntitiestoRead = SentinelString.size();
	
	TPZStack <int> GeneralData(NumEntitiestoRead,0);
	TPZStack <int> DataToProcess(NumEntitiestoRead,-1);
	
	{		
		// reading a general map information by filter
		std::ifstream read (FileName.c_str());
		std::string FlagString;
		long cont = 0;
		
		while(read)
		{
			char buf[1024];
			read.getline(buf, 1024);
			std::string str(buf);
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			if(SentinelString[cont] == "" || SentinelString[cont] == "\r") 
			{
				cont++;	
			}			
			
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if( (str != "" || str != "\r") && FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) {
					char buftemp[1024];
					read.getline(buftemp, 1024);
					std::string strtemp(buftemp);
					GeneralData[cont]++;
					if(strtemp == "" || strtemp == "\r") 
					{
						FlagString = "";
						GeneralData[cont]--;
						cont++;
						break;
					}
				}
				
			}	
			
			
		}
	}
	
	
	for (long i = 0 ; i < NumEntitiestoRead; i++ )
	{

		if(SentinelString[i] == "KABSOLUTE" || SentinelString[i] == "KABSOLUTE\r")
		{
			numKData=GeneralData[i];
			DataToProcess[i]=0;		
		}
	}
	
	numelements=numKData;
	TPZFMatrix<REAL> Kabsolute(3,3);
	Kabsolute.Resize(3,3);
	Kabsolute.Zero();
//	TPZStack<TPZFMatrix<REAL> > KabsoluteMap(numelements,0);
	TPZStack<TPZFMatrix<REAL> > KabsoluteMap(100,0);	
	
	long elementId = 0;
	long ContOfKs = 0;
	
	REAL kxx , kxy, kxz;
	REAL kyx , kyy, kyz;
	REAL kzx , kzy, kzz;	
	
	
	
	{
		
		// reading a general mesh information by filter
		std::ifstream read (FileName.c_str());
		std::string FlagString;
		long cont = 0;
		int dim = 0;
		int flag = 0;
		while(read)
		{
			char buf[1024];	
			read.getline(buf, 1024);
			std::string str(buf);
			std::string strtemp="InitialState";			
			
			// Reading General Data
			if(str == SentinelString[cont]) 
			{
				FlagString = str;	
			}
			
			if(SentinelString[cont] == "" || SentinelString[cont] == "\r") 
			{
				cont++;	
			}			
			
			if(SentinelString[cont] == "EndReading") 
			{
				break;	
			}			
			
			if( (str != "" || str != "\r" )&& FlagString == SentinelString[cont]) 
			{
				// Data scaning
				while (read) 
				{
					
					switch (DataToProcess[cont]) 
					{
						case 0:
						{
							//"KABSOLUTE"
							if (GeneralData[cont] != 0)
							{
								read >> elementId;
								read >> kxx;
								read >> kxy;
								read >> kxz;
								read >> kyx;
								read >> kyy;
								read >> kyz;
								read >> kzx;
								read >> kzy;
								read >> kzz;								
								
								Kabsolute(0,0)=kxx;
								Kabsolute(0,1)=kxy;
								Kabsolute(0,2)=kxz;
								Kabsolute(1,0)=kyx;
								Kabsolute(1,1)=kyy;
								Kabsolute(1,2)=kyz;
								Kabsolute(2,0)=kzx;								
								Kabsolute(2,1)=kzy;								
								Kabsolute(2,2)=kzz;
								
								KabsoluteMap[elementId]=Kabsolute;
								
								ContOfKs++;					
							}
							if(ContOfKs == numKData)
							{
								strtemp = "";	
							}						
						}
							break;					
						default:
						{
							strtemp = "";
						}
							break;
					}		
					
					if(strtemp == "" || strtemp == "\r") 
					{
						FlagString = "";
						cont++;
						break;
					}
				}
				
			}	
			
			
		}
	}	
	
	
	SetKMap(KabsoluteMap);
	
}



// Contribute methods

void TPZMultiphase::Contribute(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef) {
	DebugStop();
}


void TPZMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef){


#ifdef DEBUG
	int nref =  datavec.size();
	if (nref != 3 ) 
	{
        std::cout << " Erro. The size of the datavec is different from 3 \n";
		DebugStop();
	}
#endif
	
	
    // Getting weight functions	
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiP =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiS =  datavec[2].phi;	

    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFMatrix<REAL> &dphiS = datavec[2].dphix;// it is a zero value since S is constant by element.
	
	// number of test functions for each state variable
    int phrQ, phrP, phrS;
    phrQ = datavec[0].fVecShapeIndex.NElements();
    phrP = phiP.Rows();	
    phrS = phiS.Rows();
	
    //	Getting and computing another required data
    REAL TimeStep = this->fDeltaT;
    REAL Theta = this->fTheta;
    int GeoID = datavec[0].gelElId;	
    TPZFMatrix<REAL> Kabsolute;
    if (fYorN)
    {

	    Kabsolute=fKabsoluteMap[GeoID];
    }
    else 
    {
	    K(Kabsolute);
    }

	

    TPZManVector<REAL,3> sol_q =datavec[0].sol[0];
    TPZManVector<REAL,3> sol_p =datavec[1].sol[0];
    TPZManVector<REAL,3> sol_s =datavec[2].sol[0];

    REAL rockporosity, oildensity, waterdensity;
    REAL drockporositydp, doildensitydp, dwaterdensitydp;	

    REAL oilviscosity, waterviscosity;
    REAL doilviscositydp, dwaterviscositydp;

    REAL bulklambda, oillambda, waterlambda;	
    REAL dbulklambdadp, doillambdadp, dwaterlambdadp;
    REAL dbulklambdads, doillambdads, dwaterlambdads;	

    // Functions computed at point x_{k} for each integration point
    int VecPos= 0;
    REAL PressureRef = 1.0e6;
    Porosity(PressureRef, rockporosity, drockporositydp);
    RhoOil(PressureRef, oildensity, doildensitydp);
    RhoWater(PressureRef, waterdensity, dwaterdensitydp);
    OilViscosity(sol_p[VecPos], oilviscosity, doilviscositydp);
    WaterViscosity(sol_p[VecPos], waterviscosity, dwaterviscositydp);	
    OilLabmda(oillambda, sol_p[VecPos], sol_s[VecPos], doillambdadp, doillambdads);
    WaterLabmda(waterlambda, sol_p[VecPos], sol_s[VecPos], dwaterlambdadp, dwaterlambdads);
    Labmda(bulklambda, sol_p[VecPos], sol_s[VecPos], dbulklambdadp, dbulklambdads);	
	

	//	////////////////////////// Jacobian Matrix ///////////////////////////////////
	//	Contribution of domain integrals for Jacobian matrix
	// n+1 time step 
	if(gState == ECurrentState)
	{	 	
		
		//	First Block (Equation One) constitutive law
		// Integrate[(viscosity/density)*dot(v,v), Omega_{e} ]	(Equation One)
		REAL OneOverLambda = 1.0/bulklambda;		
		for(int iq=0; iq<phrQ; iq++)
		{
			
			int ivectorindex		= datavec[0].fVecShapeIndex[iq].first;
			int ishapeindex			= datavec[0].fVecShapeIndex[iq].second;			
			for (int jq=0; jq<phrQ; jq++)
			{

				int jvectorindex	= datavec[0].fVecShapeIndex[jq].first;
				int jshapeindex		= datavec[0].fVecShapeIndex[jq].second;
				
				REAL dotprod = 
				(phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex)) * (phiQ(jshapeindex,0)*datavec[0].fNormalVec(0,jvectorindex)) +
				(phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex)) * (phiQ(jshapeindex,0)*datavec[0].fNormalVec(1,jvectorindex)) +
				(phiQ(ishapeindex,0)*datavec[0].fNormalVec(2,ivectorindex)) * (phiQ(jshapeindex,0)*datavec[0].fNormalVec(2,jvectorindex)) ;	//	dot(q,v)
				
				ek(iq,jq) += weight * OneOverLambda * dotprod;
				
			}					

		}
		

		//	First Block (Equation One) constitutive law
		// Integrate[(d(1/bulklambdal)/dS)*dot(q,v), Omega_{e} ]	(Equation One)
		/*REAL OneOverLambda = 1/bulklambda;*/		
		for(int iq=0; iq<phrQ; iq++)
		{
			
			int ivectorindex		= datavec[0].fVecShapeIndex[iq].first;
			int ishapeindex			= datavec[0].fVecShapeIndex[iq].second;			
			for (int jsat=0; jsat<phrS; jsat++)
			{
				
			REAL dotprod = 
			(sol_q[0]) * (phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex)) +
			(sol_q[1]) * (phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex)) ;			
				
				ek(iq,jsat + phrQ + phrP) -= weight * dbulklambdads  * OneOverLambda * OneOverLambda * phiS(jsat,0) * dotprod;
				
			}					

		}
		
		
			
		
		for(int iq=0; iq<phrQ; iq++)
		{
			
			int ivectorindex		= datavec[0].fVecShapeIndex[iq].first;
			int ishapeindex			= datavec[0].fVecShapeIndex[iq].second;
					
			for (int jp=0; jp<phrP; jp++)
			{
				
				//	Compute grad(W)
				TPZManVector<STATE> dsolp(2,0);		
				dsolp[0] = dphiP(0,jp)*datavec[1].axes(0,0)+dphiP(1,jp)*datavec[1].axes(1,0);
				dsolp[1] = dphiP(0,jp)*datavec[1].axes(0,1)+dphiP(1,jp)*datavec[1].axes(1,1);				
				
				REAL e1e1	=	(Kabsolute(0,0)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex))+	
								 Kabsolute(0,1)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex))+
								 Kabsolute(0,2)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(2,ivectorindex)))*(dsolp[0]);
				
				REAL e2e2	=	(Kabsolute(1,0)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex))+	
								 Kabsolute(1,1)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex))+
								 Kabsolute(1,2)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(2,ivectorindex)))*(dsolp[1]);
				
					ek(iq, phrQ + jp) += weight * ( e1e1 + e2e2 );

			}
			
			
		}
		
			//	Second Block (Equation Two) Bulk flux  equation
			// Integrate[dot(grad(w),v), Omega_{e}]	(Equation Two)
			for(int ip=0; ip<phrP; ip++)
			{
				//	Compute grad(W)
				TPZManVector<STATE> dsolp(2,0);		
				dsolp[0] = dphiP(0,ip)*datavec[1].axes(0,0)+dphiP(1,ip)*datavec[1].axes(1,0);
				dsolp[1] = dphiP(0,ip)*datavec[1].axes(0,1)+dphiP(1,ip)*datavec[1].axes(1,1);			
				
				for (int jq=0; jq<phrQ; jq++)
				{
					
					int jvectorindex	= datavec[0].fVecShapeIndex[jq].first;
					int jshapeindex		= datavec[0].fVecShapeIndex[jq].second;
					
					REAL dotprod = 
					(dsolp[0]) * (phiQ(jshapeindex,0)*datavec[0].fNormalVec(0,jvectorindex)) +
					(dsolp[1]) * (phiQ(jshapeindex,0)*datavec[0].fNormalVec(1,jvectorindex)) ;
					
					ek(ip+phrQ,jq) += weight * dotprod;						
					
				}
				
			}
		
		//	Third Vector Block (Equation three)
		//  Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
		for(int isat=0; isat<phrS; isat++)
		{

			for (int jsat=0; jsat<phrS; jsat++)
			{
				ek(isat+phrQ+phrP,jsat+phrQ+phrP) +=  (rockporosity*waterdensity) * weight * phiS(isat,0) * phiS(jsat,0);
			}

		}		

		
		//	////////////////////////// Jacobian Matrix ///////////////////////////////////
		//	End of contribution of domain integrals for Jacobian matrix

	}
	

	//	n time step
	//	This values are constant in Newton iteration
	if(gState == ELastState)
	{
		
	}
	
	
	
	// Compute Residual contribution ef for domain integrals
	this->Contribute(datavec, weight, ef);

    
}

//  Residual vector contribution
void TPZMultiphase::Contribute(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef)
{
    // Getting weight functions
    TPZFMatrix<REAL>  &phiQ =  datavec[0].phi;
    TPZFMatrix<REAL>  &phiP =  datavec[1].phi;
    TPZFMatrix<REAL>  &phiS =  datavec[2].phi;	

    TPZFMatrix<REAL> &dphiQ = datavec[0].dphix;
    TPZFMatrix<REAL> &dphiP = datavec[1].dphix;
    TPZFMatrix<REAL> &dphiS = datavec[2].dphix;// This a null value since S is constant by element.

    // number of test functions for each state variable
    int phrQ, phrP, phrS;
    phrQ = datavec[0].fVecShapeIndex.NElements(); //	phiQ.Rows();
    phrP = phiP.Rows();	
    phrS = phiS.Rows();

    //	Getting and computing another required data
    int GeoID = datavec[0].gelElId;	
    TPZFMatrix<REAL> Kabsolute;
    if (fYorN) 
    {
	    Kabsolute=fKabsoluteMap[GeoID];
    }
    else 
    {
	    K(Kabsolute);
    }

	
    TPZManVector<REAL,3> sol_q =datavec[0].sol[0];
    TPZManVector<REAL,3> sol_p =datavec[1].sol[0];
    TPZManVector<REAL,3> sol_s =datavec[2].sol[0];
    TPZFMatrix<REAL> dsol_q =datavec[0].dsol[0];
    TPZFMatrix<REAL> dsol_p =datavec[1].dsol[0];
    TPZFMatrix<REAL> dsol_s =datavec[2].dsol[0];
	
    TPZFMatrix<> axesQ, axesP;
    axesQ=datavec[0].axes;
    
    REAL Pressure = sol_p[0];
    
    REAL rockporosity, oildensity, waterdensity;
    REAL drockporositydp, doildensitydp, dwaterdensitydp;	
    
    REAL oilviscosity, waterviscosity;
    REAL doilviscositydp, dwaterviscositydp;
    
    REAL bulklambda, oillambda, waterlambda;	
    REAL dbulklambdadp, doillambdadp, dwaterlambdadp;
    REAL dbulklambdads, doillambdads, dwaterlambdads;	
    
    // Functions computed at point x_{k} for each integration point
    int VecPos= 0;
    REAL PressureRef = 1.0e6;
    Porosity(PressureRef, rockporosity, drockporositydp);
    RhoOil(PressureRef, oildensity, doildensitydp);
    RhoWater(PressureRef, waterdensity, dwaterdensitydp);
    OilViscosity(sol_p[VecPos], oilviscosity, doilviscositydp);
    WaterViscosity(sol_p[VecPos], waterviscosity, dwaterviscositydp);	
    OilLabmda(oillambda, sol_p[VecPos], sol_s[VecPos], doillambdadp, doillambdads);
    WaterLabmda(waterlambda, sol_p[VecPos], sol_s[VecPos], dwaterlambdadp, dwaterlambdads);
    Labmda(bulklambda, sol_p[VecPos], sol_s[VecPos], dbulklambdadp, dbulklambdads);		
    
    //	////////////////////////// Residual Vector ///////////////////////////////////
    //	Contribution of domain integrals for Residual Vector

	// n+1 time step 
	if(gState == ECurrentState)
	{	
//		This block was verified	
		
		//	First Block (Equation One) constitutive law
		// Integrate[(viscosity/density)*dot(q,v), Omega_{e}]	(Equation One)
		REAL ViscOverdensity = waterviscosity/waterdensity;
		REAL OneOverLambda = 1.0/bulklambda;			
		for(int iq=0; iq<phrQ; iq++)
		{
			
			int ivectorindex		= datavec[0].fVecShapeIndex[iq].first;
			int ishapeindex			= datavec[0].fVecShapeIndex[iq].second;			
				
			REAL dotprod = 
			(sol_q[0]) * (phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex)) +
			(sol_q[1]) * (phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex)) ;			
			
			ef(iq) +=  OneOverLambda * weight * dotprod;
			
		}
		
		
			//		This block was verified						
			//	First Block (Equation One) constitutive law
			//	Integrate [ K dot(v,grad(P)) , Omega_{e}]	(Equation One)
			
			//	Compute grad(P)
			TPZManVector<STATE> dsolp(2,0);		
			dsolp[0] = dsol_p(0,0)*datavec[1].axes(0,0)+dsol_p(1,0)*datavec[1].axes(1,0);
			dsolp[1] = dsol_p(0,0)*datavec[1].axes(0,1)+dsol_p(1,0)*datavec[1].axes(1,1);		
			
			for(int iq=0; iq<phrQ; iq++)
			{
				
				int ivectorindex		= datavec[0].fVecShapeIndex[iq].first;
				int ishapeindex			= datavec[0].fVecShapeIndex[iq].second;					
				
				
				REAL e1e1	=	(Kabsolute(0,0)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex))+	
								 Kabsolute(0,1)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex))+
								 Kabsolute(0,2)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(2,ivectorindex)))*(dsolp[0]);
				
				REAL e2e2	=	(Kabsolute(1,0)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(0,ivectorindex))+	
								 Kabsolute(1,1)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(1,ivectorindex))+
								 Kabsolute(1,2)*(phiQ(ishapeindex,0)*datavec[0].fNormalVec(2,ivectorindex)))*(dsolp[1]);
				
				ef(iq) += weight * (e1e1 + e2e2);
				
			}		
			
			
			//	Second Block (Equation Two) Bulk flux  equation
			// Integrate[dot(grad(w),q), Omega_{e}]	(Equation Two)
	//		This block was verified		
			for(int ip=0; ip<phrP; ip++)
			{

				//	Compute grad(W)
				TPZManVector<STATE> dsolp(2,0);		
				dsolp[0] = dphiP(0,ip)*datavec[1].axes(0,0)+dphiP(1,ip)*datavec[1].axes(1,0);
				dsolp[1] = dphiP(0,ip)*datavec[1].axes(0,1)+dphiP(1,ip)*datavec[1].axes(1,1);	
				
					REAL dotprod = 
					(dsolp[0]) * (sol_q[0]) +
					(dsolp[1]) * (sol_q[1]);
				
					ef(ip+phrQ) += weight * dotprod;				
				
			}
	
		
//		This block was verified
		
		//	Third Vector Block (Equation three)
		// Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
		REAL SaturationAtnplusOne = sol_s[0];		
		for(int isat=0; isat<phrS; isat++)
		{
			ef(isat+phrQ+phrP) += weight * (rockporosity * waterdensity) * phiS(isat,0) * SaturationAtnplusOne;				
		}	
		
	}
	

	// n time step 
	if(gState == ELastState)
	{
		  
	  
//		This block was verified	
		
		//	Third Vector Block (Equation three)
		// (-1.0) * Integrate[porosity*density*L*S, Omega_{e}] (Equation three)
		REAL SaturationAtnTimeStep = sol_s[0]; //	Gettin Saturation at n time step
		for(int isat=0; isat<phrS; isat++)
		{
				ef(isat+phrQ+phrP) -= weight * (rockporosity*waterdensity) * phiS(isat,0) * SaturationAtnTimeStep;				
		}		
		

	}	
	
	//	////////////////////////// Residual Vector ///////////////////////////////////
	//	End of contribution of domain integrals for Residual Vector	
	

}


void TPZMultiphase::ContributeInterface(TPZVec<TPZMaterialData> &datavec,TPZVec<TPZMaterialData> &dataleftvec,TPZVec<TPZMaterialData> &datarightvec,REAL weight,TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMultiphase::ContributeInterface(TPZMaterialData &data, TPZMaterialData &dataleft, TPZMaterialData &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	DebugStop();
}

void TPZMultiphase::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef)
{
	
      TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
      TPZFMatrix<REAL> &phiQR = dataright[0].phi;
      TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
      TPZFMatrix<REAL> &phiPR = dataright[1].phi;	
      TPZFMatrix<REAL> &phiSL = dataleft[2].phi;
      TPZFMatrix<REAL> &phiSR = dataright[2].phi;

      TPZManVector<REAL,3> &normal = data.normal;	

      REAL n1 = normal[0];
      REAL n2 = normal[1];	
      //	REAL n3 = normal[2];

      TPZManVector<REAL,3> sol_qL =dataleft[0].sol[0];
      TPZManVector<REAL,3> sol_pL =dataleft[1].sol[0];
      TPZManVector<REAL,3> sol_sL =dataleft[2].sol[0];

      TPZManVector<REAL,3> sol_qR =dataright[0].sol[0];
      TPZManVector<REAL,3> sol_pR =dataright[1].sol[0];
      TPZManVector<REAL,3> sol_sR =dataright[2].sol[0];	


      //	Getting Q solution for left and right side
      REAL qxL = sol_qL[0];
      REAL qyL = sol_qL[1];
      REAL qxR = sol_qR[0];
      REAL qyR = sol_qR[1];	
      //	REAL qz = sol_qR[2];		

      //	Getting S solution for left and right side	
      REAL SaturationL	=	sol_sL[0];
      REAL SaturationR	=	sol_sR[0];

      //	Getting another required data
      REAL TimeStep = this->fDeltaT;
      REAL Theta = this->fTheta;

      // Getting Harmonic mean of permeabilities	
      TPZFMatrix<REAL> KabsoluteLeft;
      TPZFMatrix<REAL> KabsoluteRight;
      TPZFMatrix<REAL> Kmean(3,3,0);

      int GeoIDLeft = dataleft[0].gelElId;
      int GeoIDRight = dataright[0].gelElId;	

      REAL rockporosityl, oildensityl, waterdensityl;
      REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;	
      REAL rockporosityr, oildensityr, waterdensityr;
      REAL drockporositydpr, doildensitydpr, dwaterdensitydpr;

      REAL oilviscosityl, waterviscosityl;
      REAL doilviscositydpl, dwaterviscositydpl;
      REAL oilviscosityr, waterviscosityr;
      REAL doilviscositydpr, dwaterviscositydpr;	

      REAL bulklambdal, oillambdal, waterlambdal;	
      REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
      REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;	
      REAL bulklambdar, oillambdar, waterlambdar;	
      REAL dbulklambdadpr, doillambdadpr, dwaterlambdadpr;
      REAL dbulklambdadsr, doillambdadsr, dwaterlambdadsr;		

       REAL bulkfoill, bulkfwaterl;	
      REAL dbulkfoildpl, dbulkfwaterdpl;
      REAL dbulkfoildsl, dbulkfwaterdsl;      

      REAL bulkfoilr, bulkfwaterr;	
      REAL dbulkfoildpr, dbulkfwaterdpr;
      REAL dbulkfoildsr, dbulkfwaterdsr;      

      // Functions computed at point x_{k} for each integration point
      int VecPos= 0;
      REAL PressureRef = 1.0e6;
      Porosity(PressureRef, rockporosityl, drockporositydpl);
      RhoOil(PressureRef, oildensityl, doildensitydpl);
      RhoWater(PressureRef, waterdensityl, dwaterdensitydpl);
      OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
      WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);	
      OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
      WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
      Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
      fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
      fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        

      
      Porosity(PressureRef, rockporosityr, drockporositydpr);
      RhoOil(PressureRef, oildensityr, doildensitydpr);
      RhoWater(PressureRef, waterdensityr, dwaterdensitydpr);
      OilViscosity(sol_pR[VecPos], oilviscosityr, doilviscositydpr);
      WaterViscosity(sol_pR[VecPos], waterviscosityr, dwaterviscositydpr);	
      OilLabmda(oillambdar, sol_pR[VecPos], sol_sR[VecPos], doillambdadpr, doillambdadsr);
      WaterLabmda(waterlambdar, sol_pR[VecPos], sol_sR[VecPos], dwaterlambdadpr, dwaterlambdadsr);
      Labmda(bulklambdar, sol_pR[VecPos], sol_sR[VecPos], dbulklambdadpr, dbulklambdadsr);	
      fOil(bulkfoilr, sol_pR[VecPos], sol_sR[VecPos], dbulkfoildpr, dbulkfoildsr);
      fWater(bulkfwaterr, sol_pR[VecPos], sol_sR[VecPos], dbulkfwaterdpr, dbulkfwaterdsr);

      if (fYorN) 
      {

	  KabsoluteLeft=fKabsoluteMap[GeoIDLeft];
	  KabsoluteRight=fKabsoluteMap[GeoIDRight];
	  
      }
      else 
      {
	  K(KabsoluteLeft);
	  K(KabsoluteRight);	
      }	
	
		
	for (int in=0; in < 3; in++) {
		for (int jn=0; jn < 3; jn++) {
			if (KabsoluteLeft(in,jn)+KabsoluteRight(in,jn)==0) {
				Kmean(in,jn)= 0.0;
			}
			else 
			{
				Kmean(in,jn)= (2.0*KabsoluteLeft(in,jn)*KabsoluteRight(in,jn))/(KabsoluteLeft(in,jn)+KabsoluteRight(in,jn));
			}
			
		}
	}
	
	

//    int nsolQL = dataleft[0].sol.NElements();
//    int nsolQR = dataright[0].sol.NElements();	
	

	int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
	int QRowsRight = dataright[0].fVecShapeIndex.NElements();	
	int PRowsleft = phiPL.Rows();
	int PRowsRight = phiPR.Rows();
	int SRowsleft = phiSL.Rows();
	int SRowsRight = phiSR.Rows();	
	
	int iRightInterfaceBlock = QRowsleft + PRowsleft + SRowsleft;
	int jRightInterfaceBlock = QRowsleft + PRowsleft + SRowsleft;
	
    
	if(gState == ECurrentState)
	{	 
						
// 		if (fnewWS) 
// 		{
			//	////////////////////////// Jacobian Matrix ///////////////////////////////////
			//	Contribution of contour integrals for Jacobian matrix		
			
			
			//	First Block (Equation One) constitutive law
			// Integrate[L dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
			for (int iq=0; iq < QRowsleft; iq++) 
			{
				
				int iLvectorindex		= dataright[0].fVecShapeIndex[iq].first;
				int iLshapeindex		= dataright[0].fVecShapeIndex[iq].second;			
				
				for (int jp=0; jp<PRowsleft; jp++)
				{					
					
					REAL e1e1	=	(Kmean(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kmean(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kmean(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
					
					REAL e2e2	=	(Kmean(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kmean(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kmean(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
					
					ek(iq, QRowsleft + jp) -= weight * (e1e1 + e2e2 ) * phiPL(jp,0);				
					
				}
				
			}
			
			//	First Block (Equation One) constitutive law
			// Integrate[L dot(K v, n), Gamme_{e}]	(Equation One) Right-Right Part		
			for (int iq=0; iq < QRowsRight; iq++) 
			{
				int iRvectorindex		= dataright[0].fVecShapeIndex[iq].first;
				int iRshapeindex		= dataright[0].fVecShapeIndex[iq].second;			
				
				for (int jp=0; jp<PRowsRight; jp++)
				{	
					
					
					REAL e1e1	=	(Kmean(0,0)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(0,iRvectorindex))+	
									 Kmean(0,1)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(1,iRvectorindex))+
									 Kmean(0,2)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(2,iRvectorindex)))*(n1);
					
					REAL e2e2	=	(Kmean(1,0)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(0,iRvectorindex))+	
									 Kmean(1,1)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(1,iRvectorindex))+
									 Kmean(1,2)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(2,iRvectorindex)))*(n2);
					
					
					ek(iq + iRightInterfaceBlock, QRowsRight + jp + jRightInterfaceBlock) +=  weight * (e1e1 + e2e2 ) * phiPR(jp,0) ;
					
				}
				
			}		
			
			
			
			//	Second Block (Equation Two) Bulk flux  equation
			// Integrate[L dot(v, n), Gamma_{e}]	(Equation Two) Left-Left Part		
			for (int ip=0; ip < PRowsleft; ip++) 
			{
				
				for (int jq=0; jq<QRowsleft; jq++)
				{
					
					int jvectorindex	= dataleft[0].fVecShapeIndex[jq].first;
					int jshapeindex		= dataleft[0].fVecShapeIndex[jq].second;
					
					REAL dotprod = 
					(n1) * (phiQL(jshapeindex,0)*dataleft[0].fNormalVec(0,jvectorindex)) +
					(n2) * (phiQL(jshapeindex,0)*dataleft[0].fNormalVec(1,jvectorindex)) ;
					
					ek(ip+QRowsleft,jq) -= weight * dotprod * phiPL(ip,0);					
					
				}
				
			}
			
			
			//	Second Block (Equation Two) Bulk flux  equation
			// Integrate[L dot(v, n), Gamma_{e}]	(Equation Two) Right-Right Part		
			for (int ip=0; ip < PRowsRight; ip++) 
			{
				
				for (int jq=0; jq<QRowsRight; jq++)
				{
					
					int jvectorindex	= dataright[0].fVecShapeIndex[jq].first;
					int jshapeindex		= dataright[0].fVecShapeIndex[jq].second;
					
					REAL dotprod = 
					(n1) * (phiQR(jshapeindex,0)*dataright[0].fNormalVec(0,jvectorindex)) +
					(n2) * (phiQR(jshapeindex,0)*dataright[0].fNormalVec(1,jvectorindex)) ;
					
					ek(ip+QRowsRight+iRightInterfaceBlock,jq+jRightInterfaceBlock) += weight * dotprod * phiPR(ip,0);					
					
				}
				
			}
// 		}
	
		
		//	Upwind scheme
		//	Third Vector Block (Equation three) Saturation  equation
		REAL dotqnL = (qxL*n1) + (qyL*n2);
		REAL dotqnR = (qxR*n1) + (qyR*n2);		
		REAL UpwindSaturation = 0.0;		

		if (dotqnL > 0.0) 
		{
			UpwindSaturation = bulkfwaterl;			
			
			//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Left-Left Part			
			for(int isat=0; isat<SRowsleft; isat++)
			{
				for(int jsat=0; jsat<SRowsleft; jsat++)
				{
					ek(isat+QRowsleft+PRowsleft,jsat+QRowsleft+PRowsleft) 
					+= weight * (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
				}
			}

			//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Right-Left Part			
			for(int isat=0; isat<SRowsRight; isat++)
			{
				
				for(int jsat=0; jsat<SRowsleft; jsat++)
				{
					ek(isat+QRowsRight+PRowsRight+iRightInterfaceBlock,jsat+QRowsleft+PRowsleft) 
					-= weight * (Theta) * (TimeStep) * phiSR(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
				}				
			}	
			
		}	
		
		else 
		
		{
			UpwindSaturation = bulkfwaterr;			
			
		//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Left-Right Part			
			for(int isat=0; isat<SRowsleft; isat++)
			{
				for(int jsat=0; jsat<SRowsRight; jsat++)
				{
					ek(isat+QRowsleft+PRowsleft,jsat+QRowsRight+PRowsRight+jRightInterfaceBlock) 
					+= weight * (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsr * phiSR(jsat,0) * dotqnL;
				}
			}
		//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Right-Right Part			
			for(int isat=0; isat<SRowsRight; isat++)
			{
				
				for(int jsat=0; jsat<SRowsRight; jsat++)
				{
					ek(isat+QRowsRight+PRowsRight+iRightInterfaceBlock,jsat+QRowsRight+PRowsRight+jRightInterfaceBlock) 
					-= weight * (Theta) * (TimeStep) * phiSR(isat,0) * dbulkfwaterdsr * phiSR(jsat,0) * dotqnL;
				}				
			}
			
		}


		//	Third Vector Block (Equation three) Saturation  equation
		//	Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]	(Equation three) Left-Left Part		
		for (int isat=0; isat < SRowsleft; isat++) {
			
			for (int jq=0; jq<QRowsleft; jq++)
			{
				int jLvectorindex		= dataleft[0].fVecShapeIndex[jq].first;
				int jLshapeindex		= dataleft[0].fVecShapeIndex[jq].second;
				
				
				REAL dotprodL = 
				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(0,jLvectorindex)) * (n1) +
				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(1,jLvectorindex)) * (n2) ;//+
				//				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(2,jLvectorindex)) * (n3) ;	//	dot(q,n)	left
				
				ek(isat+QRowsleft+PRowsleft,jq)	+= weight * (Theta) * (TimeStep) * phiSL(isat,0) * UpwindSaturation * dotprodL;
				
			}
		}
		
		//	Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]	(Equation three) Right-Left Part		
		for (int isat=0; isat < SRowsRight; isat++) 
		{
			for (int jq=0; jq<QRowsleft; jq++)
			{
				int jLvectorindex		= dataleft[0].fVecShapeIndex[jq].first;
				int jLshapeindex		= dataleft[0].fVecShapeIndex[jq].second;
				
				
				REAL dotprodL = 
				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(0,jLvectorindex)) * (n1) +
				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(1,jLvectorindex)) * (n2) ;//+
				//				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(2,jLvectorindex)) * (n3) ;	//	dot(q,n)	left
				
				ek(isat+QRowsRight+PRowsRight+iRightInterfaceBlock,jq)	-= weight * (Theta) * (TimeStep) * phiSR(isat,0) * UpwindSaturation * dotprodL;
				
			}				
		}			
			
		
	}
	
	
	//	////////////////////////// Jacobian Matrix ///////////////////////////////////
	//	End of contribution of countour integrals for Jacobian matrix	
	
	
	
	//	n time step
	//	This values are constant in Newton iteration
	if(gState == ELastState)
	{
	
	}		
	
	if (fnewWS) {
		// Compute Residual contribution ef for contour integrals	
		this->ContributeInterface(data, dataleft, dataright, weight, ef);
	}


}



void TPZMultiphase::ContributeInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, TPZVec<TPZMaterialData> &dataright, REAL weight, TPZFMatrix<STATE> &ef)
{
	
      TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
      TPZFMatrix<REAL> &phiQR = dataright[0].phi;
      TPZFMatrix<REAL> &phiPL = dataleft[1].phi;
      TPZFMatrix<REAL> &phiPR = dataright[1].phi;	
      TPZFMatrix<REAL> &phiSL = dataleft[2].phi;
      TPZFMatrix<REAL> &phiSR = dataright[2].phi;

      TPZManVector<REAL,3> &normal = data.normal;


      REAL n1 = normal[0];
      REAL n2 = normal[1];	


      TPZManVector<REAL,3> sol_qL =dataleft[0].sol[0];
      TPZManVector<REAL,3> sol_pL =dataleft[1].sol[0];
      TPZManVector<REAL,3> sol_sL =dataleft[2].sol[0];

      TPZManVector<REAL,3> sol_qR =dataright[0].sol[0];
      TPZManVector<REAL,3> sol_pR =dataright[1].sol[0];
      TPZManVector<REAL,3> sol_sR =dataright[2].sol[0];	


      //	Getting Q solution for left and right side
      REAL qxL = sol_qL[0];
      REAL qyL = sol_qL[1];
      REAL qxR = sol_qR[0];
      REAL qyR = sol_qR[1];	
      //	REAL qz = sol_qR[2];		

      //	Getting P solution for left and right side	
      REAL PseudoPressureL	=	sol_pL[0];
      REAL PseudoPressureR	=	sol_pR[0];	

      //	Getting S solution for left and right side	
      REAL SaturationL	=	sol_sL[0];
      REAL SaturationR	=	sol_sR[0];

      REAL dotqnL = (qxL*n1) + (qyL*n2);
      REAL dotqnR = (qxR*n1) + (qyR*n2);


      //	Getting another required data
      REAL TimeStep = this->fDeltaT;
      REAL Theta = this->fTheta;

      // Getting Harmonic mean of permeabilities	
      TPZFMatrix<REAL> KabsoluteLeft;
      TPZFMatrix<REAL> KabsoluteRight;
      TPZFMatrix<REAL> Kmean(3,3,0);	

      int GeoIDLeft = dataleft[0].gelElId;	
      int GeoIDRight = dataright[0].gelElId;	

      REAL rockporosityl, oildensityl, waterdensityl;
      REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;	
      REAL rockporosityr, oildensityr, waterdensityr;
      REAL drockporositydpr, doildensitydpr, dwaterdensitydpr;

      REAL oilviscosityl, waterviscosityl;
      REAL doilviscositydpl, dwaterviscositydpl;
      REAL oilviscosityr, waterviscosityr;
      REAL doilviscositydpr, dwaterviscositydpr;	

      REAL bulklambdal, oillambdal, waterlambdal;	
      REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
      REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;	
      REAL bulklambdar, oillambdar, waterlambdar;	
      REAL dbulklambdadpr, doillambdadpr, dwaterlambdadpr;
      REAL dbulklambdadsr, doillambdadsr, dwaterlambdadsr;
      
      
      REAL bulkfoill, bulkfwaterl;	
      REAL dbulkfoildpl, dbulkfwaterdpl;
      REAL dbulkfoildsl, dbulkfwaterdsl;      

      REAL bulkfoilr, bulkfwaterr;	
      REAL dbulkfoildpr, dbulkfwaterdpr;
      REAL dbulkfoildsr, dbulkfwaterdsr;      

      // Functions computed at point x_{k} for each integration point
      int VecPos= 0;
      REAL PressureRef = 1.0e6;
      Porosity(PressureRef, rockporosityl, drockporositydpl);
      RhoOil(PressureRef, oildensityl, doildensitydpl);
      RhoWater(PressureRef, waterdensityl, dwaterdensitydpl);
      OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
      WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);	
      OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
      WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
      Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
      fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
      fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl);        

      
      Porosity(PressureRef, rockporosityr, drockporositydpr);
      RhoOil(PressureRef, oildensityr, doildensitydpr);
      RhoWater(PressureRef, waterdensityr, dwaterdensitydpr);
      OilViscosity(sol_pR[VecPos], oilviscosityr, doilviscositydpr);
      WaterViscosity(sol_pR[VecPos], waterviscosityr, dwaterviscositydpr);	
      OilLabmda(oillambdar, sol_pR[VecPos], sol_sR[VecPos], doillambdadpr, doillambdadsr);
      WaterLabmda(waterlambdar, sol_pR[VecPos], sol_sR[VecPos], dwaterlambdadpr, dwaterlambdadsr);
      Labmda(bulklambdar, sol_pR[VecPos], sol_sR[VecPos], dbulklambdadpr, dbulklambdadsr);	
      fOil(bulkfoilr, sol_pR[VecPos], sol_sR[VecPos], dbulkfoildpr, dbulkfoildsr);
      fWater(bulkfwaterr, sol_pR[VecPos], sol_sR[VecPos], dbulkfwaterdpr, dbulkfwaterdsr);

      if (fYorN) 
      {

	KabsoluteLeft=fKabsoluteMap[GeoIDLeft];

	KabsoluteRight=fKabsoluteMap[GeoIDRight];
      }
      else 
      {
	K(KabsoluteLeft);
	K(KabsoluteRight);
	K(Kmean);		
      }		

	
	for (int in=0; in < 3; in++) {
		for (int jn=0; jn < 3; jn++) {
			if (KabsoluteLeft(in,jn)+KabsoluteRight(in,jn)<=0) {
				Kmean(in,jn)= 0.0;
			}
			else 
			{
				Kmean(in,jn)= (2.0*KabsoluteLeft(in,jn)*KabsoluteRight(in,jn))/(KabsoluteLeft(in,jn)+KabsoluteRight(in,jn));
			}

		}
	}
	
	
	int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
	int QRowsRight = dataright[0].fVecShapeIndex.NElements();	
	int PRowsleft = phiPL.Rows();
	int PRowsRight = phiPR.Rows();
	int SRowsleft = phiSL.Rows();
	int SRowsRight = phiSR.Rows();	
	
	int iRightInterfaceBlock = QRowsleft + PRowsleft + SRowsleft;
	int jRightInterfaceBlock = QRowsleft + PRowsleft + SRowsleft;	
	
	
	//	////////////////////////// Residual Vector ///////////////////////////////////
	//	Contribution of contour integrals for Residual Vector
	
	// n+1 time step 
	if(gState == ECurrentState)
	{	
		
					
		//This block was verified			
		//	First Block (Equation One) constitutive law
		// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
		for (int iq=0; iq < QRowsleft; iq++) 
		{
			
				int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
				int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
				
				
				REAL e1e1	=	(Kmean(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
								 Kmean(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
								 Kmean(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
				
				REAL e2e2	=	(Kmean(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
								 Kmean(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
								 Kmean(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
				
				ef(iq) -= weight * (e1e1 + e2e2 ) * PseudoPressureL;
			
		}
		
		//This block was verified			
		//	First Block (Equation One) constitutive law
		// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Right-Right Part		
		for (int iq=0; iq < QRowsRight; iq++) 
		{
			
				int iRvectorindex		= dataright[0].fVecShapeIndex[iq].first;
				int iRshapeindex		= dataright[0].fVecShapeIndex[iq].second;
				
				
				REAL e1e1	=	(Kmean(0,0)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(0,iRvectorindex))+	
								 Kmean(0,1)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(1,iRvectorindex))+
								 Kmean(0,2)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(2,iRvectorindex)))*(n1);
				
				REAL e2e2	=	(Kmean(1,0)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(0,iRvectorindex))+	
								 Kmean(1,1)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(1,iRvectorindex))+
								 Kmean(1,2)*(phiQR(iRshapeindex,0)*dataright[0].fNormalVec(2,iRvectorindex)))*(n2);
		
				
				ef(iq+iRightInterfaceBlock) += weight * (e1e1 + e2e2 ) * PseudoPressureR;
			
		}	
		
		
		//	Second Block (Equation Two) Bulk flux  equation
		// Integrate[L dot(q, n), Gamme_{e}]	(Equation Two) Left-Left Part		
		for (int ip=0; ip < PRowsleft; ip++) 
		{
				ef(ip+QRowsleft) -= weight * dotqnL * phiPL(ip,0);						
		}
		
		
		//	Second Block (Equation Two) Bulk flux  equation
		// Integrate[L dot(q, n), Gamme_{e}]	(Equation Two) Right-Right Part		
		for (int ip=0; ip < PRowsRight; ip++) 
		{
				ef(ip+QRowsRight+iRightInterfaceBlock) += weight * dotqnR * phiPR(ip,0);	
		}
		
	
		REAL UpwindSaturation = 0.0;
		
		if (dotqnL > 0.0) 
		{
			UpwindSaturation = bulkfwaterl;				
			
		}
		else 
		{
			UpwindSaturation = bulkfwaterr;
			
		}			
		
		//	Third Vector Block (Equation three)
		// (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Left-Left Part
		for(int isat=0; isat<SRowsleft; isat++)
		{
			REAL ResidualPart	=	(Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
			ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
		}
		
		// (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
		for(int isat=0; isat<SRowsRight; isat++)
		{
			REAL ResidualPart	=	(Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
			ef(isat+QRowsRight+PRowsRight+iRightInterfaceBlock) -= weight * phiSR(isat,0) * ResidualPart;
		}
		
		
	}
	
	
	
	
	// n time step 
	if(gState == ELastState)
	{
		
		REAL dotqnL = (qxL*n1) + (qyL*n2);
		REAL dotqnR = (qxR*n1) + (qyR*n2);
		REAL UpwindSaturation = 0.0;
		
		if (dotqnL > 0.0) 
		{
			UpwindSaturation = bulkfwaterl;
						
		}
		else 
		{
			UpwindSaturation = bulkfwaterr;
			
		}
		
		//	Third Vector Block (Equation three)
		// (1-Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Left-Left Part
		for(int isat=0; isat<SRowsleft; isat++)
		{
			
			REAL ResidualPart	=	(1-Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
			ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
		}
		
		// (1-Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
		for(int isat=0; isat<SRowsRight; isat++)
		{
			
			REAL ResidualPart	=	(1-Theta) * (TimeStep) * ( UpwindSaturation * dotqnL );
			ef(isat+QRowsRight+PRowsRight+iRightInterfaceBlock) -= weight * phiSR(isat,0) * ResidualPart;
		}					

		
	
		
	}	
	
	//	////////////////////////// Residual Vector ///////////////////////////////////
	//	End of contribution of contour integrals for Residual Vector		
	
}



void TPZMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}

void TPZMultiphase::ContributeBC(TPZMaterialData &data, REAL weight, TPZFMatrix<STATE> &ek, TPZFMatrix<STATE> &ef, TPZBndCond &bc)
{
	DebugStop();
}	


void TPZMultiphase::ContributeBC(TPZVec<TPZMaterialData> &datavec, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc){
    
	
#ifdef DEBUG
    int nref =  datavec.size();
	if (nref != 3 ) {
        std::cout << " Erro.!! datavec tem que ser de tamanho 3 \n";
		DebugStop();
	}
#endif	


    
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZMaterialData &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	DebugStop();
}

void TPZMultiphase::ContributeBCInterface(TPZMaterialData &data, TPZVec<TPZMaterialData> &dataleft, REAL weight, TPZFMatrix<STATE> &ek,TPZFMatrix<STATE> &ef,TPZBndCond &bc)
{
	
      int nref =  dataleft.size();
      if (nref != 3) {
	      std::cout << " Error:: datavec size must to be equal to 3 \n";
	      DebugStop();
      }
      if (bc.Val2().Rows() != 4) 
      {
	      std::cout << " Error:: This material need boundary conditions for qx, qy, p (pore pressure) and s (Saturation) .\n";
	      std::cout << " give me one matrix with this form Val2(3,1).\n";
	      DebugStop();
      }

      if (bc.Val1().Rows() != 4) 
      {
	      std::cout << " Error:: This material need boundary conditions for qx, qy, p (pore pressure) and s (Saturation) .\n";
	      DebugStop();
      }	

      TPZFMatrix<REAL> &phiQL = dataleft[0].phi;
      TPZFMatrix<REAL> &phiPL = dataleft[1].phi;	
      TPZFMatrix<REAL> &phiSL = dataleft[2].phi;
      int QRowsleft = dataleft[0].fVecShapeIndex.NElements();
      int QRowsleft1 = phiQL.Rows();
      int PRowsleft = phiPL.Rows();
      int SRowsleft = phiSL.Rows();

      TPZManVector<REAL,3> &normal = data.normal;
      REAL n1 = normal[0];
      REAL n2 = normal[1];

      TPZManVector<REAL,3> sol_qL =dataleft[0].sol[0];
      TPZManVector<REAL,3> sol_pL =dataleft[1].sol[0];
      TPZManVector<REAL,3> sol_sL =dataleft[2].sol[0];


      //	Getting Q solution for left and right side
      REAL qxL = sol_qL[0];
      REAL qyL = sol_qL[1];
      REAL dotqnL = (qxL*n1) + (qyL*n2);


      //	Getting P solution for left side	
      REAL PseudoPressureL	=	sol_pL[0];

      //	Getting S solution for left side	
      REAL SaturationL	=	sol_sL[0];

      //	Getting another required data

      REAL TimeStep = this->fDeltaT;
      REAL Theta = this->fTheta;

      // Getting Harmonic mean of permeabilities	
      int GeoIDLeft = dataleft[0].gelElId;	
      TPZFMatrix<REAL> Kabsolute;	
      if (fYorN) 
      {
	      Kabsolute=fKabsoluteMap[GeoIDLeft];
      }
      else 
      {
	      K(Kabsolute);
      }	

      REAL rockporosityl, oildensityl, waterdensityl;
      REAL drockporositydpl, doildensitydpl, dwaterdensitydpl;	

      REAL oilviscosityl, waterviscosityl;
      REAL doilviscositydpl, dwaterviscositydpl;
      ;	

      REAL bulklambdal, oillambdal, waterlambdal;	
      REAL dbulklambdadpl, doillambdadpl, dwaterlambdadpl;
      REAL dbulklambdadsl, doillambdadsl, dwaterlambdadsl;	
	      

      REAL bulkfoill, bulkfwaterl;	
      REAL dbulkfoildpl, dbulkfwaterdpl;
      REAL dbulkfoildsl, dbulkfwaterdsl;      
     

      // Functions computed at point x_{k} for each integration point
      int VecPos= 0;
      REAL PressureRef = 1.0e6;
      Porosity(PressureRef, rockporosityl, drockporositydpl);
      RhoOil(PressureRef, oildensityl, doildensitydpl);
      RhoWater(PressureRef, waterdensityl, dwaterdensitydpl);
      OilViscosity(sol_pL[VecPos], oilviscosityl, doilviscositydpl);
      WaterViscosity(sol_pL[VecPos], waterviscosityl, dwaterviscositydpl);	
      OilLabmda(oillambdal, sol_pL[VecPos], sol_sL[VecPos], doillambdadpl, doillambdadsl);
      WaterLabmda(waterlambdal, sol_pL[VecPos], sol_sL[VecPos], dwaterlambdadpl, dwaterlambdadsl);
      Labmda(bulklambdal, sol_pL[VecPos], sol_sL[VecPos], dbulklambdadpl, dbulklambdadsl);
      fOil(bulkfoill, sol_pL[VecPos], sol_sL[VecPos], dbulkfoildpl, dbulkfoildsl);
      fWater(bulkfwaterl, sol_pL[VecPos], sol_sL[VecPos], dbulkfwaterdpl, dbulkfwaterdsl); 
		
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	
//	Regular Controur integrals
	
	// n+1 time step 
	if(gState == ECurrentState)
	{	

	
		//	First Block (Equation One) constitutive law
		// Integrate[L dot(K v, n), Gamma_{e}]	(Equation One) Left-Left part		
		for (int iq=0; iq < QRowsleft; iq++) 
		{
			
			for (int jp=0; jp<PRowsleft; jp++)
			{	
				int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
				int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
				
				
				REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
								 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
								 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
				
				REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
								 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
								 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
				
				ek(iq, QRowsleft + jp) -=  weight * (e1e1 + e2e2 ) * phiPL(jp,0);		
				
			}
			
		}
		
		//	Second Block (Equation Two) Bulk flux  equation
		// Integrate[L dot(v, n), Gamme_{e}]	(Equation Two) Left-Left Part		
		for (int ip=0; ip < PRowsleft; ip++) 
		{
			
			for (int jq=0; jq<QRowsleft; jq++)
			{
				
				int jvectorindex	= dataleft[0].fVecShapeIndex[jq].first;
				int jshapeindex		= dataleft[0].fVecShapeIndex[jq].second;
				
				REAL dotprod = 
				(n1) * (phiQL(jshapeindex,0)*dataleft[0].fNormalVec(0,jvectorindex)) +
				(n2) * (phiQL(jshapeindex,0)*dataleft[0].fNormalVec(1,jvectorindex)) ;
				
				ek(ip+QRowsleft,jq) -= weight * dotprod * phiPL(ip,0);
				
			}
			
		}
		
		
		//		This block was verified			
		//	First Block (Equation One) constitutive law
		// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
		for (int iq=0; iq < QRowsleft; iq++) 
		{
			
			int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
			int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
			
			
			REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
							 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
							 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
			
			REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
							 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
							 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
			
			ef(iq) -= weight * (e1e1 + e2e2) * PseudoPressureL;
			
		}
		
		//	Second Block (Equation Two) Bulk flux  equation
		// Integrate[L dot(q, n), Gamme_{e}]	(Equation Two) Left-Left Part		
		for (int ip=0; ip < PRowsleft; ip++) 
		{
			ef(ip+QRowsleft) -=  weight * dotqnL * phiPL(ip,0);						
		}	


	
	}	
			
	
	


	// n time step 
	if(gState == ELastState)
	{
		
	}	
	
	
//	Regular Controur integrals
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////	

	
	STATE v2[3];
	v2[0] = bc.Val2()(0,0);	//	qx
	v2[1] = bc.Val2()(1,0);	//	qy
	v2[2] = bc.Val2()(2,0);	//	Pressure
	v2[3] = bc.Val2()(3,0);	//	Saturation
	REAL qN = (v2[0]*n1 + v2[1]*n2);	// Normal Flux	
	

	switch (bc.Type()) {
		case 0 :	// Inflow bc: S = 1 and qn = (+ or -)g 
			
		{

		}
			break;
			
		case 1 :			// Inflow bc, P = Pin and S = Sin 
		{
			if(gState == ECurrentState)
			{
					
					//		This block was verified			
					//	First Block (Equation One) constitutive law
					// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
					for (int iq=0; iq < QRowsleft; iq++) 
					{
						
						int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
						int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;

						REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
										 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
										 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
						
						REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
										 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
										 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
						
						ef(iq) += weight * (e1e1 + e2e2 ) * (v2[2]);
						
					}
					

			//	Upwind scheme
			//	Third Vector Block (Equation three) Saturation  equation
			REAL UpwindSaturation = 0.0;
				
				if (dotqnL >= 0.0) 
				{
					UpwindSaturation = bulkfwaterl;
					
					//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Bc-Left Part			
					for(int isat=0; isat<SRowsleft; isat++)
					{
						for(int jsat=0; jsat<SRowsleft; jsat++)
						{
							ek(isat+QRowsleft+PRowsleft,jsat+QRowsleft+PRowsleft) += weight * 
							(Theta) * (TimeStep) * phiSL(isat,0)* dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
						}
					}						

				}
				
				else 
				{
					fWater(bulkfwaterl, sol_pL[VecPos], v2[3], dbulkfwaterdpl, dbulkfwaterdsl); 				  
					UpwindSaturation = bulkfwaterl;	
					
				}
								

				//	Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]	(Equation three) Right-Left Part		
				for (int isat=0; isat < SRowsleft; isat++) {
					
					for (int jq=0; jq<QRowsleft; jq++)
					{
						int jLvectorindex		= dataleft[0].fVecShapeIndex[jq].first;
						int jLshapeindex		= dataleft[0].fVecShapeIndex[jq].second;

						REAL dotprodL = 
						(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(0,jLvectorindex)) * (n1) +
						(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(1,jLvectorindex)) * (n2) ;//+
						//				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(2,jLvectorindex)) * (n3) ;	//	dot(q,n)	left
						
						ek(isat+QRowsleft+PRowsleft,jq)	+= weight * (Theta) * (TimeStep) * phiSL(isat,0) * (UpwindSaturation) * dotprodL;
						
					}
					
				}
				
				// (1-Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
				for(int isat=0; isat<SRowsleft; isat++)
				{
					REAL ResidualPart	=	(Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
					ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
				}	
				
			}			
			
			if(gState == ELastState)
			{
			
				REAL dotqnL = (qxL*n1) + (qyL*n2);
				REAL UpwindSaturation = 0.0;
				
				
				if (dotqnL > 0.0) 
				{
					UpwindSaturation = bulkfwaterl;
				}
				
				else 
				{
					fWater(bulkfwaterl, sol_pL[VecPos], v2[3], dbulkfwaterdpl, dbulkfwaterdsl); 				  
					UpwindSaturation = bulkfwaterl;
				  
					
				}
				
				
				// (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
				for(int isat=0; isat<SRowsleft; isat++)
				{
					
					REAL ResidualPart	=	(1-Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
					ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
				}				
				
				
				
			}
			
			
			
			
		}
			break;
			
		case 2 :			// Inflow bc, P = Pout and Sout
		{
			if(gState == ECurrentState)
			{
				
				//		This block was verified			
				//	First Block (Equation One) constitutive law
				// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
				for (int iq=0; iq < QRowsleft; iq++) 
				{
					
					int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
					int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
					
					
					REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
					
					REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
					
					ef(iq) += weight * (e1e1 + e2e2 ) * (v2[2]);
					
				}
				
				
				
				//	Upwind scheme
				//	Third Vector Block (Equation three) Saturation  equation
				REAL UpwindSaturation = 0.0;					

				if (dotqnL >= 0.0) 
				{
										
				UpwindSaturation = bulkfwaterl;
									
					//	Theta * TimeStep * Integrate[L L^{upwind} dot(v, n), Gamme_{e}]	(Equation three) Bc-Left Part			
					for(int isat=0; isat<SRowsleft; isat++)
					{
						for(int jsat=0; jsat<SRowsleft; jsat++)
						{
							ek(isat+QRowsleft+PRowsleft,jsat+QRowsleft+PRowsleft) += weight 
							* (Theta) * (TimeStep) * phiSL(isat,0) * dbulkfwaterdsl * phiSL(jsat,0) * dotqnL;
						}
					}						
					
				}	
				
				else 
					
				{					
					
				UpwindSaturation = 0.0;//bulkfwaterl;					
					
				if (dotqnL < 0.){ std::cout << "Boundary condition error: inflow detected in outflow boundary condition: dotqnL = " << dotqnL << "\n";}					
					
				}			
						

				//	Theta * TimeStep * Integrate[L S dot(v, n), Gamme_{e}]	(Equation three) Right-Left Part		
				for (int isat=0; isat < SRowsleft; isat++) {
					
					for (int jq=0; jq<QRowsleft; jq++)
					{
						int jLvectorindex		= dataleft[0].fVecShapeIndex[jq].first;
						int jLshapeindex		= dataleft[0].fVecShapeIndex[jq].second;
						
						
						REAL dotprodL = 
						(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(0,jLvectorindex)) * (n1) +
						(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(1,jLvectorindex)) * (n2) ;//+
						//				(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(2,jLvectorindex)) * (n3) ;	//	dot(q,n)	left
						
						ek(isat+QRowsleft+PRowsleft,jq)	+= weight * (Theta) * (TimeStep) * phiSL(isat,0) * (UpwindSaturation) * dotprodL;
						
					}
					
				}
				
				// (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
				for(int isat=0; isat<SRowsleft; isat++)
				{
					
					REAL ResidualPart	=	(Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
					ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
				}				
				
				
				
				
			}			
			
			if(gState == ELastState)
			{
				
				REAL dotqnL = (qxL*n1) + (qyL*n2);
				REAL UpwindSaturation = 0.0;
				
				
				if (dotqnL >= 0.0) 
				{
					UpwindSaturation = bulkfwaterl;
				}
				
				else 
				{	
					UpwindSaturation = 0.0;//bulkfwaterl;					
						
					if (dotqnL < 0.) {std::cout << "Boundary condition error: inflow detected in outflow boundary condition: dotqnL = " << dotqnL << "\n";					
						
					}				  
					
				}
				
				
				// (Theta) * deltat * Integrate[L*S dot(q,n), Gamma_{e}] (Equation three) Right-Left Part
				for(int isat=0; isat<SRowsleft; isat++)
				{
					
					REAL ResidualPart	=	(1-Theta) * (TimeStep) * ( (UpwindSaturation) * dotqnL );
					ef(isat+QRowsleft+PRowsleft) += weight * phiSL(isat,0) * ResidualPart;
				}				
								
				
			}	
		}	
			break;			
			
		case 3 :			// Impervious bc qn = 0 (Geological flow barrier) (P and S hold free)
		{
			if(gState == ECurrentState)
			{	
				
				//	First Block (Equation One) constitutive law
				// Integrate[L dot(K v, n), Gamma_{e}]	(Equation One) Left-Left part		
				for (int iq=0; iq < QRowsleft; iq++) 
				{
					
					for (int jp=0; jp<PRowsleft; jp++)
					{	
						int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
						int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
						
						
						REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
										 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
										 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
						
						REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
										 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
										 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
						
						ek(iq, QRowsleft + jp) += weight * (e1e1 + e2e2 ) * phiPL(jp,0);		
						
					}
					
				}

	
				//		This block was verified			
				//	First Block (Equation One) constitutive law
				// Integrate[P dot(K v, n), Gamme_{e}]	(Equation One) Left-Left part		
				for (int iq=0; iq < QRowsleft; iq++) 
				{
					
					int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
					int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
					
					
					REAL e1e1	=	(Kabsolute(0,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kabsolute(0,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kabsolute(0,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n1);
					
					REAL e2e2	=	(Kabsolute(1,0)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex))+	
									 Kabsolute(1,1)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex))+
									 Kabsolute(1,2)*(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(2,iLvectorindex)))*(n2);
					
					ef(iq) += weight * (e1e1 + e2e2) * PseudoPressureL;
					
				}				
				
				
				//	Phil's Hint
				for(int iq=0; iq < QRowsleft; iq++)
				{				
					int iLvectorindex		= dataleft[0].fVecShapeIndex[iq].first;
					int iLshapeindex		= dataleft[0].fVecShapeIndex[iq].second;
					
					REAL vni	=	(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(0,iLvectorindex)*n1)+(phiQL(iLshapeindex,0)*dataleft[0].fNormalVec(1,iLvectorindex)*n2);						
					
//					ef(iq) += weight * ( (gBigNumber * ( dotqnL - qN ) * vni ) );
//					ef(iq) += weight * ( (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN  ) * vni ) );
					ef(iq) += weight * ( (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * vni ) );					
					
					for (int jq=0; jq < QRowsleft; jq++) 
					{
						int jLvectorindex		= dataleft[0].fVecShapeIndex[jq].first;
						int jLshapeindex		= dataleft[0].fVecShapeIndex[jq].second;	
						
						REAL vnj	=	(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(0,jLvectorindex)*n1)+(phiQL(jLshapeindex,0)*dataleft[0].fNormalVec(1,jLvectorindex)*n2);							
//						ek(iq,jq) += weight * ( (gBigNumber * ( vnj ) * vni ) );
//						ek(iq,jq) += weight * ( 2.0 * (gBigNumber * ( dotqnL - qN ) * ( vnj ) * vni ) );												
						ek(iq,jq) += weight * ( 4.0 * (gBigNumber * ( dotqnL - qN ) * ( dotqnL - qN ) * ( dotqnL - qN ) * ( vnj ) * vni ) );						
					}
				}

	
			}	
			
			if(gState == ELastState)
			{	
			}
			
			
		}
			break;
	}
	
#ifdef LOG4CXX
	if(logdata->isDebugEnabled())
	{
		std::stringstream sout;
		ek.Print("ek = ",sout,EMathematicaInput);
		ef.Print("ef = ",sout,EMathematicaInput);
		LOGPZ_DEBUG(logdata,sout.str())
	}
#endif	
	

}



/** Returns the variable index associated with the name */
int TPZMultiphase::VariableIndex(const std::string &name){
	if(!strcmp("BulkVelocity",name.c_str()))        return  1;
	if(!strcmp("WaterVelocity",name.c_str()))        return  2;
	if(!strcmp("OilVelocity",name.c_str()))        return  3;
	if(!strcmp("WeightedPressure",name.c_str()))    return  4;
	if(!strcmp("WaterSaturation",name.c_str()))    return  5;
	if(!strcmp("OilSaturation",name.c_str()))    return  6;
	if(!strcmp("WaterDensity",name.c_str()))    return  7;
	if(!strcmp("OilDensity",name.c_str()))    return  8;
	if(!strcmp("RockPorosity",name.c_str()))    return  9;
	
	return TPZMaterial::VariableIndex(name);
}

int TPZMultiphase::NSolutionVariables(int var){
	if(var == 1) return 2;
	if(var == 2) return 2;
	if(var == 3) return 2;
	if(var == 4) return 1;
	if(var == 5) return 1;
	if(var == 6) return 1;
	if(var == 7) return 1;
	if(var == 8) return 1;
	if(var == 9) return 1;
	
	return TPZMaterial::NSolutionVariables(var);
}

void TPZMultiphase::Solution(TPZVec<TPZMaterialData> &datavec, int var, TPZVec<STATE> &Solout){
	
	Solout.Resize( this->NSolutionVariables(var));
	TPZVec<STATE> SolP, SolQ, SolS;
	SolQ = datavec[0].sol[0];
	SolP = datavec[1].sol[0];
	SolS = datavec[2].sol[0];
	
	if(var == 1){ //function (state variable Q)
		Solout[0] = SolQ[0];
        Solout[1] = SolQ[1];
		return;
	}
    
    if(var == 4){
		Solout[0] = SolP[0];//function (state variable p)
		return;
	}
    
	if(var == 5){
		Solout[0] = SolS[0];//function (state variable S)
		return;
	}

	if(var == 6){
		Solout[0] = 1.0-SolS[0];//function (state variable S)
		return;
	}	

}


void TPZMultiphase::FillDataRequirements(TPZVec<TPZMaterialData > &datavec)
{
	int nref = datavec.size();
	for(int i = 0; i<nref; i++ )
	{
		datavec[i].SetAllRequirements(false);
		datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNeighborSol = true;
		datavec[i].fNeedsNeighborCenter = false;
		datavec[i].fNeedsNormal = true;
	}
	
}

void TPZMultiphase::FillBoundaryConditionDataRequirement(int type,TPZVec<TPZMaterialData > &datavec){
    int nref = datavec.size();
	for(int i = 0; i<nref; i++)
	{
        datavec[i].fNeedsSol = true;
		datavec[i].fNeedsNormal = true;
		datavec[i].fNeedsNeighborSol = true;		
	}
}