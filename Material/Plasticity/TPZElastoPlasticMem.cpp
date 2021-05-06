//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "TPZElastoPlasticMem.h"
#ifdef FIX_PLASTIC_TRANSLATORS
#include "TPZElastoPlasticMemTranslator.h"
#endif
#include "pzadmchunk.h"

TPZElastoPlasticMem::TPZElastoPlasticMem(): m_sigma(), m_elastoplastic_state(), m_plastic_steps(0),m_phi(0.), m_u(3,0.)
{
}

TPZElastoPlasticMem::TPZElastoPlasticMem(const TPZElastoPlasticMem & other):
m_sigma(other.m_sigma), m_elastoplastic_state(other.m_elastoplastic_state), m_plastic_steps(other.m_plastic_steps), m_phi(other.m_phi), m_u(other.m_u),m_ER(other.m_ER) {
    
}


TPZElastoPlasticMem::~TPZElastoPlasticMem(){
    
}

void TPZElastoPlasticMem::Write(TPZStream &buf, int withclassid) const
{
    m_sigma.Write(buf, withclassid);
    m_elastoplastic_state.Write(buf, withclassid);
    buf.Write(&m_plastic_steps);
    buf.Write(m_u);
    m_ER.Write(buf, withclassid);
}

void TPZElastoPlasticMem::Read(TPZStream &buf, void *context)
{
    m_sigma.Read(buf, context);
    m_elastoplastic_state.Read(buf, context);
    buf.Read(&m_plastic_steps);
    buf.Read(m_u);
    m_ER.Read(buf, context);
}

void TPZElastoPlasticMem::Print(std::ostream &out)const
{
    out << Name();
    out << "\nm_sigma = " << m_sigma;
    out << "\nm_elastoplastic_state = " << m_elastoplastic_state;
    out << "\nm_plastic_steps = " << m_plastic_steps;
    out << "\nm_u = " << m_u;
    out << "\nm_phi = " << m_phi;
    m_ER.Print(out);
}

const std::string TPZElastoPlasticMem::Name()const
{
    return "TPZElastoPlasticMem";
}

int TPZElastoPlasticMem::ClassId() const{
    return Hash("TPZElastoPlasticMem");
}

const TPZElastoPlasticMem & TPZElastoPlasticMem::operator=(const TPZElastoPlasticMem & other)
{
    
    /// check for self-assignment
    if(&other == this){
        return *this;
    }
    
    m_sigma = other.m_sigma;
    m_elastoplastic_state = other.m_elastoplastic_state;
    m_plastic_steps = other.m_plastic_steps;
    m_phi  = other.m_phi;
    m_u = other.m_u;
    m_ER = other.m_ER;
    
    return *this;
}

#ifdef FIX_PLASTIC_TRANSLATORS
template class TPZRestoreClassWithTranslator<TPZElastoPlasticMem, TPZElastoPlasticMemTranslator>;
#endif
template class TPZRestoreClass<TPZAdmChunkVector<TPZElastoPlasticMem>>;

