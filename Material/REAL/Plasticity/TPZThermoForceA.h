// $Id: TPZThermoForceA.h,v 1.9 2010-05-24 02:56:18 erick Exp $

#ifndef TPZTHERMOFORCEA_H
#define TPZTHERMOFORCEA_H

#include "pzvec.h"
#include "pzreal.h"
#include "TPZStream.h"
#include "TPZSavable.h"

/**
Classe que implementa o calculo da forca termodinamica (Souza Neto p. 144) e suas derivadas
Neste caso utiliza-se encruamento linear
*/
class TPZThermoForceA : public TPZSavable {
public:

    TPZThermoForceA() : fSigmaYield0(0), fK(0)
    {
    }
    
    public:
int ClassId() const override;


    const char * Name() const
    {
	   return "TPZThermoForceA";	
    }
	
    void Print(std::ostream & out) const
    {
       out << Name();
	   out << "\n fSigmaYield0 = " << fSigmaYield0;
	   out << "\n fK = " << fK;
    }
	
    void SetUp(REAL yield, REAL k) 
    {
       fSigmaYield0 = yield;
       fK           = k;
    }

    /**
    Calculo do valor da forca termo dinamica
    */
    template <class T>
    T Compute(const T & alpha) const;

    /**
    Calculo da derivada da forca termodinamica
    */
    template <class T>
    T ComputeTangent(const T & alpha) const;

    void Write(TPZStream &buf, int withclassid) const override{
        buf.Write(&fSigmaYield0);
        buf.Write(&fK);
    }

    void Read(TPZStream& buf, void* context) override {
        buf.Read(&fSigmaYield0);
        buf.Read(&fK);
    }

public:

    /**
    Tensao de plastificao inicial
    */
    REAL fSigmaYield0;
    /**
    Valor proporcional de encruamento
    */
    REAL fK;
};

template < class T >
T TPZThermoForceA::Compute(const T & alpha) const
{
    //coesao * exp^(-alpha/kc)^(p)
  //  T val = T(fSigmaYield0)*(exp(alpha/T(0.1))*exp(alpha/T(0.1)));
  //  T val2 =  val;
  //  return val2;
//	T var = T(fK) * alpha;
//    return T(fSigmaYield0) + var;
     return T(fSigmaYield0) + T(fK) * alpha;
}

template < class T >
T TPZThermoForceA::ComputeTangent(const T & alpha) const
{
   return fK;
}

#endif //TPZTHERMOFORCEA_H
