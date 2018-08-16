//$Id: pzporoelastoplasticmem.cpp,v 1.3 2009-10-06 01:00:06 erick Exp $

#include "pzporoelastoplasticmem.h"


TPZPoroElastoPlasticMem::TPZPoroElastoPlasticMem(): TPZElastoPlasticMem(), 
						fPorePressure(0), fdPorePressure(3,0.) { }
	
TPZPoroElastoPlasticMem::TPZPoroElastoPlasticMem(const TPZPoroElastoPlasticMem & source):
                         TPZElastoPlasticMem(source), fPorePressure(source.fPorePressure),
						 fdPorePressure(source.fdPorePressure)
{
}


TPZPoroElastoPlasticMem::~TPZPoroElastoPlasticMem(){ }

void TPZPoroElastoPlasticMem::Write(TPZStream &buf, int withclassid) const
{
	TPZElastoPlasticMem::Write(buf,withclassid);

	buf.Write(&fPorePressure,1);
	buf.Write(&fdPorePressure[0],3);
    buf.Write(&fv[0],3);
}

void TPZPoroElastoPlasticMem::Read(TPZStream &buf, void *context)
{
	TPZElastoPlasticMem::Read(buf,context);
	buf.Read(&fPorePressure,1);
	fdPorePressure.Resize(3);
	buf.Read(&fdPorePressure[0],3);
    fv.Resize(3);
    buf.Read(&fv[0],3);
}

void TPZPoroElastoPlasticMem::Print(std::ostream &out)const
{
	out << Name();
	out << "\n fPorePressure = " << fPorePressure;
	out << "\n fdPorePressure = " << fdPorePressure;
    out << "\n fv = " << fv;
	out << "\n Parent Class Data:";
	TPZElastoPlasticMem::Print(out);
}

const std::string TPZPoroElastoPlasticMem::Name()const
{
	return "TPZPoroElastoPlasticMem";	
}

int TPZPoroElastoPlasticMem::ClassId() const{
    return Hash("TPZPoroElastoPlasticMem") ^ TPZElastoPlasticMem::ClassId() << 1;
}

const TPZPoroElastoPlasticMem & TPZPoroElastoPlasticMem::operator=(const TPZPoroElastoPlasticMem & source)
{
	TPZElastoPlasticMem::operator=(source);
	fPorePressure = source.fPorePressure;
	fdPorePressure = source.fdPorePressure;
	
	return *this;
}
