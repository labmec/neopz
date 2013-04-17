//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "pzelastoplasticmem.h"
#include "pzmaterialid.h"
#include "poroelastoplasticid.h"


TPZElastoPlasticMem::TPZElastoPlasticMem(): fSigma(), fPlasticState(), fPlasticSteps(0),fPhi(0.), fDisplacement(3,0.) { }
	
TPZElastoPlasticMem::TPZElastoPlasticMem(const TPZElastoPlasticMem & source):
                          fSigma(source.fSigma), fPlasticState(source.fPlasticState), fPlasticSteps(source.fPlasticSteps), fPhi(source.fPhi), fDisplacement(source.fDisplacement) { }


TPZElastoPlasticMem::~TPZElastoPlasticMem(){ }

void TPZElastoPlasticMem::Write(TPZStream &buf, int withclassid)
{
	buf.Write(&fSigma[0],6);
	
	buf.Write(&fPlasticState.fEpsT[0],6);
	buf.Write(&fPlasticState.fEpsP[0],6);
	buf.Write(&fPlasticState.fAlpha,1);
	buf.Write(&fPlasticSteps,1);
    buf.Write(&fDisplacement[0],3);
//	buf.Write(&fPlasticState.fPhi,1);
}

void TPZElastoPlasticMem::Read(TPZStream &buf, void *context)
{
	buf.Read(&fSigma[0],6);
	
	buf.Read(&fPlasticState.fEpsT[0],6);
	buf.Read(&fPlasticState.fEpsP[0],6);
	buf.Read(&fPlasticState.fAlpha,1);
	buf.Read(&fPlasticSteps,1);
    buf.Read(&fDisplacement[0],3);
//	buf.Read(&fPlasticState.fPhi,1);
}

void TPZElastoPlasticMem::Print(std::ostream &out)const
{
	out << Name();
	out << "\nfSigma = " << fSigma;
	out << "\nfPlasticState = " << fPlasticState;
	out << "\nfPlasticSteps = " << fPlasticSteps;
    out << "\nfDisplacement = " << fDisplacement;
	out << "\nfPhi = " << fPhi;
}

const std::string TPZElastoPlasticMem::Name()const
{
	return "TPZElastoPlasticMem";	
}

const int TPZElastoPlasticMem::ClassId()const
{
	return TPZELASTOPLASTICMEM_ID;
}

const TPZElastoPlasticMem & TPZElastoPlasticMem::operator=(const TPZElastoPlasticMem & source)
{
	fSigma = source.fSigma;
	fPlasticState = source.fPlasticState;
	fPlasticSteps = source.fPlasticSteps;
	fPhi  = source.fPhi;
    fDisplacement = source.fDisplacement;
	
	return *this;
}
