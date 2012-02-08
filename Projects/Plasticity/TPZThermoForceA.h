// $Id: TPZThermoForceA.h,v 1.9 2010-05-24 02:56:18 erick Exp $

#ifndef TPZTHERMOFORCEA_H
#define TPZTHERMOFORCEA_H

#include "pzvec.h"
#include "pzreal.h"

/**
Classe que implementa o calculo da forca termodinamica (Souza Neto p. 144) e suas derivadas
Neste caso utiliza-se encruamento linear
*/
class TPZThermoForceA {
public:

    TPZThermoForceA() : fSigmaYield0(0), fK(0)
    {
    }

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
	T var = T(fK) * alpha; 
   return T(fSigmaYield0) + var;
}

template < class T >
T TPZThermoForceA::ComputeTangent(const T & alpha) const
{
   return fK;
}

#endif //TPZTHERMOFORCEA_H
