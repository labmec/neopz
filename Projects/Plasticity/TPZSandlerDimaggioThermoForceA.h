// $Id: TPZSandlerDimaggioThermoForceA.h,v 1.5 2009-06-29 22:54:01 erick Exp $

#ifndef TPZSANDLERDIMAGGIOTHERMOFORCEA_H
#define TPZSANDLERDIMAGGIOTHERMOFORCEA_H

#include "pzvec.h"
#include "pzreal.h"

/**
Classe que implementa o calculo da forca termodinamica (Souza Neto p. 144) e suas derivadas
Esta implementacao considera que A=alpha. Embora esteja conceitualmente
errado, a implementacao do criterio de plastificacao com A=alpha
conduz a um modelo correto.
*/
class TPZSandlerDimaggioThermoForceA {
public:

    TPZSandlerDimaggioThermoForceA() 
    {
    }
	
    TPZSandlerDimaggioThermoForceA(const TPZSandlerDimaggioThermoForceA & source)
    {
    }

    TPZSandlerDimaggioThermoForceA & operator=(const TPZSandlerDimaggioThermoForceA & source)
    {
		return *this;
    }

	const char * Name() const
    {
	   return "TPZSandlerDimaggioThermoForceA";	
    }
	
    void Print(std::ostream & out) const
    {
		out << "\n" << this->Name();
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

private:

};

template < class T >
T TPZSandlerDimaggioThermoForceA::Compute(const T & alpha) const
{
   return alpha; 
}

template < class T >
T TPZSandlerDimaggioThermoForceA::ComputeTangent(const T & alpha) const
{
   return 1.;
}

#endif //TPZSANDLERDIMAGGIOTHERMOFORCEA_H
