//$Id: pzelastoplasticmem.cpp,v 1.6 2009-06-22 00:55:14 erick Exp $

#include "TPZElastoPlasticMem.h"
#include "TPZElastoPlasticMemTranslator.h"


TPZElastoPlasticMem::TPZElastoPlasticMem(): m_sigma(), m_elastoplastic_state(), m_plastic_steps(0),m_phi(0.), m_u(3,0.) { }

TPZElastoPlasticMem::TPZElastoPlasticMem(const TPZElastoPlasticMem & source):
m_sigma(source.m_sigma), m_elastoplastic_state(source.m_elastoplastic_state), m_plastic_steps(source.m_plastic_steps), m_phi(source.m_phi), m_u(source.m_u) { }


TPZElastoPlasticMem::~TPZElastoPlasticMem(){
    
}

void TPZElastoPlasticMem::Write(TPZStream &buf, int withclassid) const
{
    m_sigma.Write(buf, withclassid);
    m_elastoplastic_state.Write(buf, withclassid);
    buf.Write(&m_plastic_steps);
    buf.Write(m_u);
    //    buf.Write(&m_elastoplastic_state.m_phi,1);
}

void TPZElastoPlasticMem::Read(TPZStream &buf, void *context)
{
    m_sigma.Read(buf, context);
    m_elastoplastic_state.Read(buf, context);
    buf.Read(&m_plastic_steps);
    buf.Read(m_u);
    //    buf.Read(&m_elastoplastic_state.m_phi,1);
}

void TPZElastoPlasticMem::Print(std::ostream &out)const
{
    out << Name();
    out << "\nm_sigma = " << m_sigma;
    out << "\nm_elastoplastic_state = " << m_elastoplastic_state;
    out << "\nm_plastic_steps = " << m_plastic_steps;
    out << "\nm_u = " << m_u;
    out << "\nm_phi = " << m_phi;
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
    m_sigma = source.m_sigma;
    m_elastoplastic_state = source.m_elastoplastic_state;
    m_plastic_steps = source.m_plastic_steps;
    m_phi  = source.m_phi;
    m_u = source.m_u;
    
    return *this;
}

//void TPZElastoPlasticMem::SetElastoPlasticState(const TPZPlasticState<REAL> & elastoplastic_state){
//    m_elastoplastic_state = elastoplastic_state;
//}
//
//
//TPZPlasticState<REAL> & TPZElastoPlasticMem::GetElastoPlasticState(){
//    return m_elastoplastic_state;
//}
//
//void TPZElastoPlasticMem::SetSigma(const TPZTensor<REAL> & sigma){
//    m_sigma = sigma;
//}
//
//TPZTensor<REAL> & TPZElastoPlasticMem::GetSigma(){
//    return m_sigma;
//}
//
//void TPZElastoPlasticMem::SetDisplacement(const TPZManVector<REAL,3> & u){
//    m_u = u;
//}
//
//TPZManVector<REAL,3> & TPZElastoPlasticMem::GetDisplacement(){
//    return m_u;
//}
//
//void TPZElastoPlasticMem::SetPlasticSteps(const int & plastic_steps){
//    m_plastic_steps = plastic_steps;
//}
//
//int & TPZElastoPlasticMem::GetPlasticSteps(){
//    return m_plastic_steps;
//}
//
//void TPZElastoPlasticMem::SetPhi(const REAL & phi){
//    m_phi = phi;
//}
//
//REAL & TPZElastoPlasticMem::GetPhi(){
//    return m_phi;
//}


template class TPZRestoreClassWithTranslator<TPZElastoPlasticMem, TPZElastoPlasticMemTranslator>;
template class TPZRestoreClass<TPZAdmChunkVector<TPZElastoPlasticMem>>;

