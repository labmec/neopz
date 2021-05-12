// $Id: TPZPlasticState.h,v 1.8 2009-07-19 05:47:18 erick Exp $

#ifndef TPZPLASTICSTATE_H
#define TPZPLASTICSTATE_H

#include "fadType.h"
#include "TPZTensor.h"
#include <iostream>

/**
 * This class holds the complete set of state variables to define a valid elastoplastic strain state
 */
template <class T>
class TPZPlasticState : public TPZSavable {
    
public:
    
    /// Tensors representing the total and plastic strain states
    TPZTensor<T> m_eps_t, m_eps_p;
    
    /// Plastic volumetric hardeing variable
    T m_hardening;
    
    /// Identifier for the regime of the material behaviour
    int m_m_type;
    
public:
    
    /// Default constructor - all values set to zero
    TPZPlasticState(): m_eps_t(), m_eps_p(), m_hardening(T(0.)), m_m_type(0) { }
    
    /// Constructor enabling predefinition of hardening
    TPZPlasticState(const T & hardening):m_eps_t(T(0.)), m_eps_p(T(0.)), m_hardening(hardening), m_m_type(0) { }
    
    /// Copy constructor
    TPZPlasticState(const TPZPlasticState<T> & source):
    m_eps_t(source.m_eps_t), m_eps_p(source.m_eps_p), m_hardening(source.m_hardening), m_m_type(source.m_m_type){ }
    
    /// Destructor
    ~TPZPlasticState(){ }
    
    /// Operator =
    const TPZPlasticState<T> & operator=(const TPZPlasticState<T> & source);
    
    /// Operator -=
    const TPZPlasticState<T> & operator-=(const TPZPlasticState<T> & source);
    
    /// Operator +=
    const TPZPlasticState<T> & operator+=(const TPZPlasticState<T> & source);
    
    /// Operator *=
    const TPZPlasticState<T> & operator*=(const TPZPlasticState<T> & source);
    
    /// Operator<<
    friend std::ostream& operator<<( std::ostream& Out, const TPZPlasticState<T> & s )
    {
        s.Print(Out);
        return Out;
    }
    
    /// More complete then Operator << because it allows derivatives supression.
    void Print(std::ostream& Out, int fadDerivatives = 1)const;
    
    // class Access Members (needed for const PlasticState access)
    
    const TPZTensor<T> & EpsT() const
    { return m_eps_t; }
    
    const TPZTensor<T> & EpsP() const
    { return m_eps_p; }
    
    const T & VolHardening() const
    { return m_hardening; }
    
    const int & MType() const
    { return m_m_type; }
        int ClassId() const override;
    
    void Read(TPZStream& buf, void* context) override {
        m_eps_t.Read(buf,context);
        m_eps_p.Read(buf,context);
        
        buf.Read(&m_hardening);
        buf.Read(&m_m_type);
    }
    
    void Write(TPZStream &buf, int withclassid) const override{
        m_eps_t.Write(buf,withclassid);
        m_eps_p.Write(buf,withclassid);
        
        buf.Write(&m_hardening);
        buf.Write(&m_m_type);
    }
    
    void CleanUp() {
        m_eps_t.Zero();
        m_eps_p.Zero();
        m_hardening = T(0.);
        m_m_type = 0;
    }
    
    /**
     * Similar to Operator=, but allows copies among different template specializations.
     * When using FAD class Types as template classes, NO DERIVATIVES are copied using this functions,
     * enabling copies from REAL to FAD types and also vice-versa.
     */
    template <class T1>
    void CopyTo(TPZPlasticState<T1> & target) const;
    
};

template <class T>
int TPZPlasticState<T>::ClassId() const{
    return Hash("TPZPlasticState") ^ ClassIdOrHash<T>() << 1;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator=(const TPZPlasticState<T> & source)
{
    m_eps_t = source.EpsT();
    m_eps_p = source.EpsP();
    m_hardening = source.VolHardening();
    m_m_type = source.MType();
    
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator+=(const TPZPlasticState<T> & source)
{
    m_eps_t += source.EpsT();
    m_eps_p += source.EpsP();
    m_hardening+= source.VolHardening();
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator-=(const TPZPlasticState<T> & source)
{
    m_eps_t -= source.EpsT();
    m_eps_p -= source.EpsP();
    m_hardening-= source.VolHardening();
    return *this;
}

template <class T>
inline const TPZPlasticState<T> & TPZPlasticState<T>::operator*=(const TPZPlasticState<T> & source)
{
    m_eps_t *= source.EpsT();
    m_eps_p *= source.EpsP();
    m_hardening*= source.VolHardening();
    return *this;
}

template <class T>
inline void TPZPlasticState<T>::Print(std::ostream& Out, int fadDerivatives)const {
    if (fadDerivatives) {
        Out << "\tm_eps_t = ";
        for (int i = 0; i < 6; i++)Out << m_eps_t[i] << " ";
        Out << std::endl;
        Out << "\tm_eps_p = ";
        for (int i = 0; i < 6; i++)Out << m_eps_p[i] << " ";
        Out << std::endl;
        Out << "\tm_hardening = " << m_hardening << std::endl;
    } else {
        Out << "\tm_eps_t = ";
        for (int i = 0; i < 6; i++)Out << TPZExtractVal::val(m_eps_t[i]) << " ";
        Out << std::endl;
        Out << "\tm_eps_p = ";
        for (int i = 0; i < 6; i++)Out << TPZExtractVal::val(m_eps_p[i]) << " ";
        Out << std::endl;
        Out << "\tm_hardening = " << TPZExtractVal::val(m_hardening) << std::endl;
    }
    Out << "\tm_m_type = " << m_m_type << std::endl;
}

template <class T>
template <class T1>
void TPZPlasticState<T>::CopyTo(TPZPlasticState<T1> & target) const
{
    EpsT().CopyTo(target.m_eps_t);
    EpsP().CopyTo(target.m_eps_p);
    target.m_hardening = TPZExtractVal::val( VolHardening() );
    target.m_m_type = MType();
}


#endif

