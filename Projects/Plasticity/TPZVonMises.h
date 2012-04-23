/*
 *  TPZVonMises.h
 *  ElastoPlasticModels
 *
 *  Created by Diogo Cecilio on 12/17/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#ifndef TPZVONMISES_H
#define TPZVONMISES_H

#include "pzlog.h"
#include "TPZPlasticStep.h"
#include "TPZYCVonMises.h"
#include "TPZThermoForceA.h"
#include "TPZElasticResponse.h"
#include "pzvec_extras.h"
#include "pzsave.h"
#include "TPZPlasticStepID.h"

#define VONMISESPARENT TPZPlasticStep<TPZYCVonMises, TPZThermoForceA, TPZElasticResponse>


class TPZVonMises : public VONMISESPARENT, public TPZSaveable  {
	
public:
	
	enum {NYield =1};
	
public:
	
    TPZVonMises():VONMISESPARENT()
    {
		fMaterialTensionSign  = 1; // internally in this material tension is negative
		fInterfaceTensionSign =  1; // by default
    }
	
    TPZVonMises(const TPZVonMises & source):VONMISESPARENT(source)
    {
    }
	
    TPZVonMises & operator=(const TPZVonMises & source)
    {
		VONMISESPARENT::operator=(source);		
		return *this;
    }
	
	virtual const char * Name() const
	{
		return "TPZVonMises";	
	}
	
	
	void SetUp()
	{
		
		//VONMISESPARENT::fYC.SetUp();
		//VONMISESPARENT::fTFA.SetUp(210./* sigmaY = MPa */,    1164./*endurecimento para cada 0.1% de deformacao em MPa(290/0.249)*/);
		//VONMISESPARENT::fER.SetUp(/*young*/ 210000., /*poisson*/ 0.3);
		
	}
	
	REAL YieldRadius(TPZPlasticState<REAL> state)
	{
				REAL radius = VONMISESPARENT::fTFA.fSigmaYield0;
				return radius;
	}
	
	
	virtual void Print(std::ostream & out) const
	{
		out << "\n" << this->Name();
		out << "\n Base Class Data:\n";
		VONMISESPARENT::Print(out);		
	}
	
	virtual int ClassId() const
	{	
		return TPZVONMISES_ID;	
	}
	
	virtual void Write(TPZStream &buf, int withclassid)
	{
	}
	
	virtual void Read(TPZStream &buf, void *context)
	{
	}	
    
    static void Steel(TPZVonMises & material)
    {
        material.fTFA.SetUp(210./* sigmaY = MPa */,    1164./*endurecimento para cada 0.1% de deformacao em MPa(290/0.249)*/);
		material.fER.SetUp(/*young*/ 210000., /*poisson*/ 0.3);
    }
	
	
public:
	
};


#endif //TPZVonMises_H
