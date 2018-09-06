//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "TPZElastoPlasticMem.h"
#include "TPZElastoPlasticMemTranslator.h"


TPZElastoPlasticMem::TPZElastoPlasticMem(): fSigma(), fPlasticState(), fPlasticSteps(0),fPhi(0.), fDisplacement(3,0.) { }
	
TPZElastoPlasticMem::TPZElastoPlasticMem(const TPZElastoPlasticMem & source):
                          fSigma(source.fSigma), fPlasticState(source.fPlasticState), fPlasticSteps(source.fPlasticSteps), fPhi(source.fPhi), fDisplacement(source.fDisplacement) { }


TPZElastoPlasticMem::~TPZElastoPlasticMem(){ }

void TPZElastoPlasticMem::Write(TPZStream &buf, int withclassid) const
{   
    fSigma.Write(buf, withclassid);
    fPlasticState.Write(buf, withclassid);
    buf.Write(&fPlasticSteps);
    buf.Write(fDisplacement);
//	buf.Write(&fPlasticState.fPhi,1);
}

void TPZElastoPlasticMem::Read(TPZStream &buf, void *context)
{
    fSigma.Read(buf, context);
    fPlasticState.Read(buf, context);
    buf.Read(&fPlasticSteps);
    buf.Read(fDisplacement);
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

int TPZElastoPlasticMem::ClassId() const{
    return Hash("TPZElastoPlasticMem");
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

template class TPZRestoreClassWithTranslator<TPZElastoPlasticMem, TPZElastoPlasticMemTranslator>;
template class TPZRestoreClass<TPZAdmChunkVector<TPZElastoPlasticMem>>;
