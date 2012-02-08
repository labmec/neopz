//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "pzelastoplasticmem.h"
#include "pzmaterialid.h"
#include "poroelastoplasticid.h"


TPZElastoPlasticMem::TPZElastoPlasticMem(): fSigma(), fPlasticState(), fPlasticSteps(0),fPhi() { }
	
TPZElastoPlasticMem::TPZElastoPlasticMem(const TPZElastoPlasticMem & source):
                          fSigma(source.fSigma), fPlasticState(source.fPlasticState), fPlasticSteps(source.fPlasticSteps), fPhi(source.fPhi) { }


TPZElastoPlasticMem::~TPZElastoPlasticMem(){ }

void TPZElastoPlasticMem::Write(TPZStream &buf, int withclassid)
{
	buf.Write(&fSigma.fData[0],6);
	
	buf.Write(&fPlasticState.fEpsT.fData[0],6);
	buf.Write(&fPlasticState.fEpsP.fData[0],6);
	buf.Write(&fPlasticState.fAlpha,1);
	buf.Write(&fPlasticSteps,1);
//	buf.Write(&fPlasticState.fPhi,1);
}

void TPZElastoPlasticMem::Read(TPZStream &buf, void *context)
{
	buf.Read(&fSigma.fData[0],6);
	
	buf.Read(&fPlasticState.fEpsT.fData[0],6);
	buf.Read(&fPlasticState.fEpsP.fData[0],6);
	buf.Read(&fPlasticState.fAlpha,1);
	buf.Read(&fPlasticSteps,1);
//	buf.Read(&fPlasticState.fPhi,1);
}

void TPZElastoPlasticMem::Print(std::ostream &out)const
{
	out << Name();
	out << "\n fSigma = " << fSigma;
	out << "\n fPlasticState = " << fPlasticState;
	out << "\n fPlasticSteps = " << fPlasticSteps;
	out << "\n fPhi = " << fPhi;
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
	
	return *this;
}
