/**
 * @file
 * @brief Contains implementations of the TPZConservationLaw methods.
 */

#include "pzconslaw.h" 
#include "pzelmat.h"
#include "pzbndcond.h"
#include "pzmatrix.h"
#include "pzfmatrix.h"
#include "pzerror.h"
#include <math.h>
#include <stdlib.h>
using namespace std;


TPZConservationLaw::TPZConservationLaw(int nummat,REAL timeStep,int dim) :
TPZDiscontinuousGalerkin(nummat),
fDim(dim),
fTimeStep(0),
fCFL(0),
fGamma(1.4),
fContributionTime(Last_CT),
fResidualType(Flux_RT)

{
	fTimeStep = timeStep;
	if(timeStep < 0 || timeStep > 1)
	{
		PZError << "TPZConservationLaw::TPZConservationLaw time step parameter > 1 , default 1.0\n";
		fTimeStep = 1.0;
	}
	
	if(dim < 1 || dim > 3)
	{
		PZError << "TPZConservationLaw::TPZConservationLaw (abort) error dimension = " << dim << endl;
		DebugStop();
	}
	fDim = dim;
	fResidualType = Residual_RT;
	
}


TPZConservationLaw::~TPZConservationLaw()
{
}

void TPZConservationLaw::Print(std::ostream &out)
{
	out << "name of material : " << Name() << "\n";
	out << "properties : \n";
	out << "\tdimension: " << fDim << endl;
	out << "\ttime step: " << fTimeStep << endl;
	out << "\tCFL: " << fCFL << endl;
	//   out << "\tDelta (diffusive term): " << fDelta << endl;
	out << "\tGamma: " << fGamma << endl;
	TPZMaterial::Print(out);
	
	
	switch(fContributionTime)
	{
		case Advanced_CT:
			out << "Advanced Contribution\n";
			break;
		case Last_CT:
			out << "Last state Contribution\n";
			break;
		case None_CT:
			out << "Contribution undefined\n";
	}
}

void TPZConservationLaw::Write(TPZStream &buf, int withclassid)
{
	TPZMaterial::Write(buf, withclassid);
	buf.Write(&fDim,1);
	buf.Write(&fTimeStep,1);
	buf.Write(&fCFL,1);
	buf.Write(&fGamma,1);
}

void TPZConservationLaw::Read(TPZStream &buf, void *context)
{
	TPZMaterial::Read(buf, context);
	buf.Read(&fDim,1);
	buf.Read(&fTimeStep,1);
	buf.Read(&fCFL,1);
	buf.Read(&fGamma,1);
	
	fContributionTime = Last_CT;
	fResidualType = Residual_RT;
}

